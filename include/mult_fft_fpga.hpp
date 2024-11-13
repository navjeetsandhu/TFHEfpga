#pragma once
#include <array>
#include "mult_fft.hpp"
#include "INTorus.hpp"
#include "params.hpp"
#include "utils.hpp"
#ifdef USE_FPGA
#include <fft_processor_fpga.h>

template <int N>
inline void TwistFpgaFFT(std::array<uint64_t, N> &res, const std::array<double, N> &a)
{
    fftFpgaLvl1.execute_direct_torus64(res.data(), a.data());
}

template <int N>
inline void TwistFpgaFFT(std::array<uint32_t, N> &res, const std::array<double, N> &a)
{
    fftFpgaLvl1.execute_direct_torus32(res.data(), a.data());
}

template <int N>
inline void TwistFpgaFFTrescale(std::array<uint64_t, N> &res, const std::array<double, N> &a, const double delta)
{
    fftFpgaLvl1.execute_direct_torus32_rescale(res.data(), a.data(), delta);
}

template <int N>
inline void TwistFpgaFFTrescale(std::array<uint32_t, N> &res, const std::array<double, N> &a, const double delta)
{
    fftFpgaLvl1.execute_direct_torus64_rescale(res.data(), a.data(), delta);
}

template <int N>
inline void TwistFpgaIFFT(std::array<double, N> &res, const std::array<uint64_t, N> &a)
{
    fftFpgaLvl1.execute_reverse_torus64(res.data(), a.data());
}

template <int N>
inline void TwistFpgaIFFT(std::array<double, N> &res, const std::array<uint32_t, N> &a)
{
    fftFpgaLvl1.execute_reverse_torus32(res.data(), a.data());
}

namespace TFHEpp {
template <class P>
inline void TwistFpgaFFTrescale(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            fftFpgaLvl1.execute_direct_torus32_rescale(res.data(), a.data(),
                                                       P::delta);
        else if constexpr (std::is_same_v<typename P::T, uint64_t>)
            fftFpgaLvl1.execute_direct_torus64_rescale(res.data(), a.data(),
                                                       P::delta);
        else
            static_assert(false_v<typename P::T>, "TwistFpgaFFTrescale!");
    }
    else
        static_assert(false_v<typename P::T>, "TwistFpgaFFTrescale!!");
}
}

#else

template <int N>
inline void TwistFpgaFFT(std::array<uint64_t, N> &res, const std::array<double, N> &a)
{
    fftplvl1.execute_direct_torus64(res.data(), a.data());
}

template <int N>
inline void TwistFpgaFFT(std::array<uint32_t, N> &res, const std::array<double, N> &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}


template <int N>
inline void TwistFpgaIFFT(std::array<double, N> &res, const std::array<uint64_t, N> &a)
{
    fftplvl1.execute_reverse_torus64(res.data(), a.data());
}

template <int N>
inline void TwistFpgaIFFT(std::array<double, N> &res, const std::array<uint32_t, N> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}

namespace TFHEpp {
template <class P>
inline void TwistFpgaFFTrescale(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>) {
        if constexpr (std::is_same_v<typename P::T, uint32_t>)
            fftplvl1.execute_direct_torus32_rescale(res.data(), a.data(),
                                                    P::delta);
        else if constexpr (std::is_same_v<typename P::T, uint64_t>)
            fftplvl1.execute_direct_torus64_rescale(res.data(), a.data(),
                                                    P::delta);
        else
            static_assert(false_v<typename P::T>, "TwistFpgaFFTrescale!");
    }
    else
        static_assert(false_v<typename P::T>, "TwistFpgaFFTrescale!!");
}
}

#endif

template <class T, int N>
void PolyMulFpgaFFT(std::array<T, N> &res, const std::array<T, N>  &a,
                const std::array<T, N>  &b)
{
    alignas(64) std::array<double, N> ffta{};
    TwistFpgaIFFT<N>(ffta, a);
    alignas(64) std::array<double, N> fftb{};
    TwistFpgaIFFT<N>(fftb, b);
    MulInFD<N>(ffta, fftb);
    TwistFpgaFFT<N>(res, ffta);
}




