#include "simulator.hpp"
#include <cmath>
#include <random>
#include <fftw3.h>
#include <algorithm>

const float G = 9.81f;
const float PI = 3.1415926535f;

OceanSimulator::OceanSimulator(int N, float L, float A, float wind_v, float wind_dir) 
    : N(N), L(L), A(A) {
    
    int size = N * N;
    // メモリ確保（高さ + 水平変位x, y）
    h0_tilde = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);
    h_tilde_t = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);
    dx_tilde_t = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);
    dy_tilde_t = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);

    height_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);
    dx_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);
    dy_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * size);

    // 3つの逆FFTプラン作成
    plan_h  = fftwf_plan_dft_2d(N, N, h_tilde_t, height_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_dx = fftwf_plan_dft_2d(N, N, dx_tilde_t, dx_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_dy = fftwf_plan_dft_2d(N, N, dy_tilde_t, dy_out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // 風の向きベクトル
    float wx = wind_v * cos(wind_dir);
    float wy = wind_v * sin(wind_dir);

    // 初期スペクトル h0 の生成
    std::default_random_engine generator;
    std::normal_distribution<float> dist(0.0, 1.0);

    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            int idx = n * N + m;

            float kn = (n <= N / 2) ? (float)n : (float)n - N;
            float km = (m <= N / 2) ? (float)m : (float)m - N;
            
            float kx = (2.0f * PI * kn) / L;
            float ky = (2.0f * PI * km) / L;
            float k_mag = sqrt(kx * kx + ky * ky);

            if (k_mag < 0.0001f) {
                h0_tilde[idx][0] = h0_tilde[idx][1] = 0.0f;
                continue;
            }

            float L_val = (wind_v * wind_v) / G;
            float phillips = A * (exp(-1.0f / (k_mag * L_val * k_mag * L_val)) / pow(k_mag, 4.0f)) 
                             * pow(std::abs((kx/k_mag)*(wx/wind_v) + (ky/k_mag)*(wy/wind_v)), 2.0f);

            float r1 = dist(generator);
            float r2 = dist(generator);
            h0_tilde[idx][0] = r1 * sqrt(phillips * 0.5f);
            h0_tilde[idx][1] = r2 * sqrt(phillips * 0.5f);
        }
    }
}

// update関数
void OceanSimulator::update(float t, float loop_period) {
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < N; ++m) {
            int idx = n * N + m;
            
            float kn = (n <= N / 2) ? (float)n : (float)n - N;
            float km = (m <= N / 2) ? (float)m : (float)m - N;
            
            float kx = (2.0f * PI * kn) / L;
            float ky = (2.0f * PI * km) / L;
            float k_mag = sqrt(kx * kx + ky * ky);
            float omega = sqrt(G * k_mag);

            // ループ周期に合わせた周波数の量子化 (Phase Snapping)
            if (loop_period > 0.0f) {
                float delta_omega = 2.0f * PI / loop_period;
                omega = std::round(omega / delta_omega) * delta_omega;
            }

            float c = cos(omega * t);
            float s = sin(omega * t);
                
            // 共役な波のインデックス (-k)
            int c_idx = ((N - n) % N) * N + ((N - m) % N);

            // --- リニアなループを実現する複素数計算 ---
            // h(k, t) = h0(k)*exp(i*wt) + conj(h0(-k))*exp(-i*wt)
            float a = h0_tilde[idx][0];
            float b = h0_tilde[idx][1];
            float c_val = h0_tilde[c_idx][0];
            float d = h0_tilde[c_idx][1];

            // Real(h) = (a+c)*cos - (b+d)*sin
            // Imag(h) = (a-c)*sin + (b-d)*cos
            float res_re = (a + c_val) * c - (b + d) * s;
            float res_im = (a - c_val) * s + (b - d) * c;

            h_tilde_t[idx][0] = res_re;
            h_tilde_t[idx][1] = res_im;

            if (k_mag < 0.0001f) {
                dx_tilde_t[idx][0] = dx_tilde_t[idx][1] = 0.0f;
                dy_tilde_t[idx][0] = dy_tilde_t[idx][1] = 0.0f;
            } else {
                float eig_x = kx / k_mag;
                float eig_y = ky / k_mag;
                // 水平変位も共役対称性を維持するように計算 (i * h)
                dx_tilde_t[idx][0] = res_im * eig_x;
                dx_tilde_t[idx][1] = -res_re * eig_x;
                dy_tilde_t[idx][0] = res_im * eig_y;
                dy_tilde_t[idx][1] = -res_re * eig_y;
            }
        }
    }

    fftwf_execute(plan_h);
    fftwf_execute(plan_dx);
    fftwf_execute(plan_dy);
}

OceanSimulator::~OceanSimulator() {
    fftwf_destroy_plan(plan_h);
    fftwf_destroy_plan(plan_dx);
    fftwf_destroy_plan(plan_dy);
    fftwf_free(h0_tilde);
    fftwf_free(h_tilde_t);
    fftwf_free(dx_tilde_t);
    fftwf_free(dy_tilde_t);
    fftwf_free(height_out);
    fftwf_free(dx_out);
    fftwf_free(dy_out);
}