#include "mult_fft.hpp"
#include "mult_fft_fpga.hpp"
#include"utils2.hpp"
#include <numeric>
#include <cmath>

template <int nbits>
void test_fft(const std::array<uint64_t, 1 << nbits>& p1)
{

    std::string string_msg = "Input p1";
    constexpr int N = 1 << nbits;

    print_results<uint64_t>(string_msg, p1.data(), p1.size());

    std::array<uint64_t, N> result{};
    std::array<uint64_t, N> result1{};
    std::fill(result.begin(), result.end(), 0);

    alignas(64) std::array<double, N> fft{};
    TwistFpgaIFFT<N>(fft, p1);

    string_msg = "TwistIFFT 64 bit";
    print_results<double>(string_msg,  fft.data(), fft.size());

    TwistFpgaFFT<N>(result, fft);
    string_msg = "TwistFFT 64 bit";
    print_results<int64_t>(string_msg,  reinterpret_cast<int64_t*>(result.data()), result.size());
}


template <int nbits>
void test_fft(const std::array<uint32_t, 1 << nbits>& p1)
{

    std::string string_msg = "Input p1";
    constexpr int N = 1 << nbits;

    print_results<uint32_t>(string_msg, p1.data(), p1.size());

    std::array<uint32_t, N> result{};
    std::fill(result.begin(), result.end(), 0);

    alignas(64) std::array<double, N> fft{};
    TwistFpgaIFFT<N>(fft, p1);


    string_msg = "TwistIFFT 32 bit";
    print_results<double>(string_msg,  fft.data(), fft.size());

    TwistFpgaFFT<N>(result, fft);
    string_msg = "TwistFFT 32 bit";
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());
}



template <class P, int nbits>
void test_fft_p()
{
    constexpr int N = 1 << nbits;

    std::array<P,  N> p1{};
    std::array<P,  N> p2{};
    std::iota(p1.begin(), p1.end(), 1);
    test_fft<nbits>(p1);
}

int main()
{
    constexpr int nbit = 10;
    test_fft_p<uint32_t,nbit>();
    test_fft_p<uint64_t,nbit>();
    return 0;
}