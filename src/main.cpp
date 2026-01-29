#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include "ocean/simulator.hpp"
#include "utils/image_writer.hpp"

// N        : グリッドの解像度（2の乗数）
// L        : シミュレーション範囲
// A        : 振幅（波の高さを決める係数）
// wind_v   : 風速
// wind_dir : 波の進行方向（ラジアン表記）

const int N = 512;
const float L = 200.0f;
const float A = 5.0f; 
const float wind_v = 50.0f;
const float wind_dir = 3.1415f / 4.0f;


int main()
{
    OceanSimulator ocean(N, L, A, wind_v, wind_dir);

    const int total_frames = 60; 
    const float dt = 0.1f; // Web側の frameDuration (100ms) に合わせる
    
    for (int frame = 0; frame < total_frames; ++frame)
    {
        float t = frame * dt;
        const float loop_period = total_frames * dt; // ループの総時間
        ocean.update(t, loop_period);

        std::vector<float> height_data(N * N), norm_r(N * N), norm_g(N * N), norm_b(N * N);
        std::vector<float> foam_data(N * N);

        for (int y = 0; y < N; ++y)
        {
            for (int x = 0; x < N; ++x)
            {
                int i = y * N + x;
                int next_x = (x + 1) % N;
                int next_y = (y + 1) % N;

                // 1. 高さの取得
                height_data[i] = ocean.height_out[i][0];

                // 2. 法線の計算 (有限差分法)
                float dx = (ocean.height_out[y * N + next_x][0] - ocean.height_out[i][0]);
                float dy = (ocean.height_out[next_y * N + x][0] - ocean.height_out[i][0]);

                // ベクトル (-dx, -dy, 1.0) を正規化
                float mag = sqrt(dx * dx + dy * dy + 1.0f);
                norm_r[i] = -dx / mag;
                norm_g[i] = -dy / mag;
                norm_b[i] = 1.0f / mag;

                // 3. ジャコビアン (Foam / 白波) の計算
                // 水平変位の微分から「メッシュがどれだけ圧縮されているか」を出す
                float dDx_dx = (ocean.dx_out[y * N + next_x][0] - ocean.dx_out[i][0]);
                float dDy_dy = (ocean.dy_out[next_y * N + x][0] - ocean.dy_out[i][0]);
                float J = (1.0f + dDx_dx) * (1.0f + dDy_dy);
                foam_data[i] = (J < 0.5f) ? 1.0f : 0.0f; // 閾値以下なら白波
            }
        }

        // テクスチャの保存
        save_heightmap("output/height_" + std::to_string(frame + 1) + ".png", N, N, height_data);
        save_image_rgb("output/normal_" + std::to_string(frame + 1) + ".png", N, N, norm_r, norm_g, norm_b);
        save_heightmap("output/foam_" + std::to_string(frame + 1) + ".png", N, N, foam_data);

        std::cout << "Saved frame " << frame + 1 << " (Height, Normal, Foam)" << std::setw(10)
                  << "(" << frame + 1 << "/60)" << std::endl;
    }

    std::cout << "Done!" << std::endl;

    return 0;
}