#include <iostream>
#include <math.h>
#ifdef USE_FPGA
#include "fpga.h"
#else
#include "fpga_test.h"
#endif

using namespace std;


/**
 * \brief  create random single precision complex floating point values
 * \param  inp : pointer to float2 data of size N
 * \param  N   : number of points in the array
 * \return true if successful
 */
void create_data(const unsigned num, float2 *inp){
    if(inp == NULL || num < 4){ throw "Bad args in create data function";}

    for(unsigned i = 0; i < num; i++){
        inp[i].x = (float)((float)rand() / (float)RAND_MAX);
        inp[i].y = (float)((float)rand() / (float)RAND_MAX);
    }
}

std::ostream& operator<<(std::ostream& os, const fpga_timing& p) {
    os << "exec: " << p.exec_t << std::endl;
    os << "pcie_read: " << p.pcie_read_t << std::endl;
    os << "pcie_write: " << p.pcie_write_t << std::endl;
    os << "svm_copyin: " << p.svm_copyin_t << std::endl;
    os << "svm_copyout: " << p.svm_copyout_t << std::endl;
    os << "valid: " << p.valid << std::endl;
    return os;
}


void print_fft_fpga_data(const unsigned num, const float2 *inp, const float2 *out)
{
    printf("\nFFT FPGA result: \n");
    for (unsigned i = 0; i < num; i++) {
        printf("%u : input - (%e %e) output - (%e %e)\n", i, inp[i].x, inp[i].y,
               out[i].x, out[i].y);
    }
    printf("\n\n");
}


void fft_data_test(const unsigned num)
{
    fpga_initialize();

    float2 *inp = new float2[num]();
    float2 *out1 = new float2[num]();
    float2 *out2 = new float2[num]();
    create_data(num, inp);

    bool inv = true;
    cout << endl << "inv: " << inv << endl;
    fpga_timing timing = fpga_fft(num, inp, out1, inv);
    cout << timing;
    print_fft_fpga_data(num, inp, out1);

    inv = false;
    cout << endl << "inv: " << inv << endl;
    timing = fpga_fft(num, out1, out2, inv);
    cout << timing;
    print_fft_fpga_data(num, out1, out2);


    delete[] inp;
    delete[] out1;
    delete[] out2;

    fpga_close();
}




int main()
{
    fft_data_test(512);
    return 0;
}