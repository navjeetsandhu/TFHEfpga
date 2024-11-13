#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>
#include "fftfpga.h"

class FFT_Processor_FPGA {
public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;

private:
    std::vector<std::complex<double>> twist;
    float2 *inbuf;
    float2 *outbuf;
    fpga_t runTimeRc;

public:
    FFT_Processor_FPGA(const int32_t N);

    void execute_reverse_int(double *res, const int32_t *a);

    void execute_reverse_torus32(double *res, const uint32_t *a);

    void execute_direct_torus32(uint32_t *res, const double *a);

    void execute_direct_torus32_rescale(uint32_t *res, const double *a,
                                        const double delta);

    void execute_reverse_torus64(double *res, const uint64_t *a);

    void execute_direct_torus64(uint64_t *res, const double *a);

    void execute_direct_torus64_rescale(uint64_t *res, const double *a,
                                        const double delta);

    ~FFT_Processor_FPGA();
};

extern FFT_Processor_FPGA fftFpgaLvl1;
