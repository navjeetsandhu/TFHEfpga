#include "mult_fft.hpp"
#include "mult_fft_fpga.hpp"
#include"utils2.hpp"
#include <numeric>
#include <cmath>

template <int nbits>
void test_mult(const std::array<uint64_t, 1 << nbits>& p1, const std::array<uint64_t, 1 << nbits>& p2)
{

    constexpr int N = 1 << nbits;
    std::string string_msg = "Naive 64 bit Multiplication";
    print_results<uint64_t>("Input p1", p1.data(), p1.size());
    print_results<uint64_t>("Input p2", p2.data(), p1.size());
    std::array<uint64_t, N> result = {};
    PolyMulNaive<uint64_t, N>(result, p1, p2);
    print_results<int64_t>(string_msg,  reinterpret_cast<int64_t*>(result.data()), result.size());


    std::fill(result.begin(), result.end(), 0);
    PolyMulFFT<uint64_t, N>(result, p1, p2);
    string_msg = "FFT 64 bit Multiplication";
    print_results<int64_t>(string_msg,  reinterpret_cast<int64_t*>(result.data()), result.size());

    std::fill(result.begin(), result.end(), 0);
    PolyMulFpgaFFT<uint64_t, N>(result, p1, p2);
    string_msg = "FPGA FFT 64 bit Multiplication";
    print_results<int64_t>(string_msg,  reinterpret_cast<int64_t*>(result.data()), result.size());
}

template <int nbits>
void test_mult(const std::array<uint32_t, 1 << nbits>& p1, const std::array<uint32_t, 1 << nbits>& p2)
{

    constexpr int N = 1 << nbits;
    std::string string_msg = "Naive 32 bit Multiplication";
    print_results<uint32_t>("Input p1", p1.data(), p1.size());
    print_results<uint32_t>("Input p2", p2.data(), p1.size());
    std::array<uint32_t, N> result = {0,0,0,0};
    PolyMulNaive<uint32_t, N>(result, p1, p2);
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());

    std::fill(result.begin(), result.end(), 0);
    PolyMulFFT<uint32_t, N>(result, p1, p2);
    string_msg = "FFT 32 bit Multiplication";
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());

    std::fill(result.begin(), result.end(), 0);
    PolyMulFpgaFFT<uint32_t, N>(result, p1, p2);
    string_msg = "FPGA FFT 32 bit Multiplication";
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());
}


template <class P, int nbits>
void test_mult()
{
    constexpr int N = 1 << nbits;

    std::array<P,  N> p1{};
    std::array<P,  N> p2{};
    std::iota(p1.begin(), p1.end(), 1);
    std::iota(p2.begin(), p2.end(), 10);
    test_mult<nbits>(p1, p2);
}

int main()
{
    constexpr int nbit = 10;
    test_mult<uint32_t,nbit>();
    test_mult<uint64_t, nbit>();
    return 0;
}