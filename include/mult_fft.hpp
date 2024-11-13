#pragma once
#include <array>
#ifdef USE_FFTW3
#include <fft_processor_fftw.h>
#else
#include <fft_processor_spqlios.h>
#endif


template <int N>
inline void TwistFFT(std::array<uint64_t, N> &res, const std::array<double, N> &a)
{
    fftplvl1.execute_direct_torus64(res.data(), a.data());
}

template <int N>
inline void TwistFFT(std::array<uint32_t, N> &res, const std::array<double, N> &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}

template <int N>
inline void TwistIFFT(std::array<double, N> &res, const std::array<uint64_t, N> &a)
{
    fftplvl1.execute_reverse_torus64(res.data(), a.data());
}

template <int N>
inline void TwistIFFT(std::array<double, N> &res, const std::array<uint32_t, N> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}

template <int N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = res[i + N / 2] * b[i + N / 2];
        double arebim = res[i] * b[i + N / 2];
        res[i] = std::fma(res[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(res[i + N / 2], b[i], arebim);
    }
}


template <class T, int N>
void PolyMulFFT(std::array<T, N> &res, const std::array<T, N>  &a,
                const std::array<T, N>  &b)
{
    alignas(64) std::array<double, N> ffta{};
    TwistIFFT<N>(ffta, a);
    alignas(64) std::array<double, N> fftb{};
    TwistIFFT<N>(fftb, b);
    MulInFD<N>(ffta, fftb);
    TwistFFT<N>(res, ffta);
}


template <class T, int N>
void PolyMulNaive(std::array<T, N> &res, const std::array<T, N>  &a,
                     const std::array<T, N>  &b){
    for (int i = 0; i < N; i++) {
        T ri = 0;
        for (int j = 0; j <= i; j++)
            ri += static_cast<typename std::make_signed<T>::type>(
                          a[j]) *
                  b[i - j];
        for (int j = i + 1; j < N; j++)
            ri -= static_cast<typename std::make_signed<T>::type>(
                          a[j]) *
                  b[N + i - j];
        res[i] = ri;
    }
}


