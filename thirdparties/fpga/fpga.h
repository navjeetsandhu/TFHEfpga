#pragma once
#include "fftfpga.h"

void fpga_initialize();
void fpga_close();
fpga_t fpga_fft(const unsigned num, const float2 *inp, float2 *out, const bool inv);
unsigned bit_reversed(unsigned x, const unsigned bits);
void correct_data_order(float2* fpgaOut, const unsigned num);