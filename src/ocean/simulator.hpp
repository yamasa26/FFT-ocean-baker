#pragma once
#include <vector>
#include <fftw3.h>

class OceanSimulator {
public:
    OceanSimulator(int N, float L, float A, float wind_v, float wind_dir);
    ~OceanSimulator();
    
    void update(float t, float loop_period);

    int N;
    float L, A;
    fftwf_complex *h0_tilde, *h_tilde_t;
    fftwf_complex *height_out;
    
    fftwf_complex *dx_tilde_t, *dy_tilde_t;
    fftwf_complex *dx_out, *dy_out;
    fftwf_plan plan_h, plan_dx, plan_dy;
};