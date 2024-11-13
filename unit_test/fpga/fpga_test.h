#pragma once
#include <iostream>
/**
 * Single Precision Complex Floating Point Data Structure
 */
typedef struct {
    float x; /**< real value */
    float y; /**< imaginary value */
} float2;

/**
 * Double Precision Complex Floating Point Data Structure
 */
typedef struct {
    double x; /**< real value */
    double y; /**< imaginary value */
} double2;


typedef struct fpga_timing {
    double pcie_read_t;     /**< Time to read from DDR to host using PCIe bus  */
    double pcie_write_t;    /**< Time to write from DDR to host using PCIe bus */
    double exec_t;          /**< Kernel execution time */
    double svm_copyin_t;    /**< Time to copy in data to SVM */
    double svm_copyout_t;   /**< Time to copy data out of SVM */
    bool valid;             /**< Represents true signifying valid execution */
} fpga_t;

fpga_t fpga_fft(const unsigned num, const float2 *inp, float2 *out, const bool inv)
{
    unsigned i;
    std::cout << std::endl << "inv: " << inv <<std::endl;
    for (i=0; i < num; i++) {
        out[i].x = inp[i].x;
        out[i].y = inp[i].y;
    }
    return {1, 2, 3, 4,5, true};
}

void fpga_initialize() {
}


void fpga_close() {
}