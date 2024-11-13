#pragma once
#include <memory>
#include "INTorus.hpp"
#include "cuhe++.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "mult_fft_fpga.hpp"
#ifdef USE_FFTW3
#include <fft_processor_fftw.h>
#else
#include <fft_processor_spqlios.h>
#endif


namespace TFHEpp {
    inline const std::unique_ptr<
        const std::array<std::array<cuHEpp::INTorus, lvl1param::n>, 2>>
            ntttwistlvl1 = cuHEpp::TwistGen<lvl1param::nbit>();
    inline const std::unique_ptr<
        const std::array<std::array<cuHEpp::INTorus, lvl1param::n>, 2>>
            ntttablelvl1 = cuHEpp::TableGen<lvl1param::nbit>();
    inline const std::unique_ptr<
            const std::array<std::array<cuHEpp::INTorus, lvl2param::n>, 2>>
            ntttwistlvl2 = cuHEpp::TwistGen<lvl2param::nbit>();
    inline const std::unique_ptr<
            const std::array<std::array<cuHEpp::INTorus, lvl2param::n>, 2>>
            ntttablelvl2 = cuHEpp::TableGen<lvl2param::nbit>();
    inline const std::unique_ptr<
            std::array<std::array<raintt::SWord, lvl1param::n>, 2>>
            raintttwist = raintt::TwistGen<lvl1param::nbit, 3>();
    inline const std::unique_ptr<
            std::array<std::array<std::array<raintt::SWord, lvl1param::n>, 2>, 2>>
            raintttable = raintt::TableGen<lvl1param::nbit>();



template <class P>
inline void TwistNTT(Polynomial<P> &res, PolynomialNTT<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
        cuHEpp::TwistNTT<typename lvl1param::T, lvl1param::nbit>(
            res, a, (*ntttablelvl1)[0], (*ntttwistlvl1)[0]);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        cuHEpp::TwistNTT<typename lvl2param::T, lvl2param::nbit>(
            res, a, (*ntttablelvl2)[0], (*ntttwistlvl2)[0]);
    else
        static_assert(false_v<typename P::T>, "Undefined TwistNTT!");
}


template <class P>
inline void TwistFFT(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
        TwistFpgaFFT<P::n>(res, a);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_direct_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistFFT!");
}


template <class P>
inline void TwistFFTrescale(Polynomial<P> &res, const PolynomialInFD<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
        TwistFpgaFFTrescale<P>(res, a);
    else if constexpr (std::is_same_v<P, lvl2param>)
        fftplvl2.execute_direct_torus64_rescale(res.data(), a.data(), P::delta);
    else
        static_assert(false_v<typename P::T>, "Undefined TwistFFT!");
}



template <class P>
inline void TwistINTT(PolynomialNTT<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
        cuHEpp::TwistINTT<typename P::T, P::nbit>(res, a, (*ntttablelvl1)[1],
                                                  (*ntttwistlvl1)[1]);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        cuHEpp::TwistINTT<typename P::T, P::nbit>(res, a, (*ntttablelvl2)[1],
                                                  (*ntttwistlvl2)[1]);
    else
        static_assert(false_v<typename P::T>, "Undefined TwistINTT!");
}

template <class P>
inline void TwistIFFT(PolynomialInFD<P> &res, const Polynomial<P> &a)
{
    if constexpr (std::is_same_v<P, lvl1param>)
        TwistFpgaIFFT<P::n>(res, a);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_reverse_torus64(res.data(), a.data());
    else
        static_assert(false_v<typename P::T>, "Undefined TwistIFFT!");
}
template <uint32_t N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        double aimbim = res[i + N / 2] * b[i + N / 2];
        double arebim = res[i] * b[i + N / 2];
        res[i] = std::fma(res[i], b[i], -aimbim);
        res[i + N / 2] = std::fma(res[i + N / 2], b[i], arebim);
    }
}

template <uint32_t N>
inline void MulInFD(std::array<double, N> &res, const std::array<double, N> &a,
                    const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = a[i] * b[i];
        res[i + N / 2] = a[i + N / 2] * b[i];
    }
    for (int i = 0; i < N / 2; i++) {
        res[i + N / 2] += a[i] * b[i + N / 2];
        res[i] -= a[i + N / 2] * b[i + N / 2];
    }
}

// Be careful about memory accesss (We assume b has relatively high memory
// access cost)
template <uint32_t N>
inline void FMAInFD(std::array<double, N> &res, const std::array<double, N> &a,
                    const std::array<double, N> &b)
{
    for (int i = 0; i < N / 2; i++) {
        res[i] = std::fma(a[i], b[i], res[i]);
        res[i + N / 2] = std::fma(a[i + N / 2], b[i], res[i + N / 2]);
    }
    for (int i = 0; i < N / 2; i++) {
        res[i + N / 2] = std::fma(a[i], b[i + N / 2], res[i + N / 2]);
        res[i] -= a[i + N / 2] * b[i + N / 2];
    }
}

template <class P>
inline void PolyMulNtt(Polynomial<P> &res, const Polynomial<P> &a,
                        const Polynomial<P> &b)
{
    PolynomialNTT<P> ntta, nttb;
    TwistINTT<P>(ntta, a);
    TwistINTT<P>(nttb, b);
    for (int i = 0; i < P::n; i++) ntta[i] *= nttb[i];
    TwistNTT<P>(res, ntta);
}
template <class P>
inline void PolyMulNaive(Polynomial<P> &res, const Polynomial<P> &a,
                         const Polynomial<P> &b)
{
    for (int i = 0; i < P::n; i++) {
        typename P::T ri = 0;
        for (int j = 0; j <= i; j++)
            ri += static_cast<typename std::make_signed<typename P::T>::type>(
                      a[j]) *
                  b[i - j];
        for (int j = i + 1; j < P::n; j++)
            ri -= static_cast<typename std::make_signed<typename P::T>::type>(
                      a[j]) *
                  b[P::n + i - j];
        res[i] = ri;
    }
}

template <class P>
inline void PolyMulFFT(Polynomial<P> &res, const Polynomial<P> &a,
                       const Polynomial<P> &b)
{
    alignas(64) PolynomialInFD<P> ffta;
    TwistIFFT<P>(ffta, a);
    alignas(64) PolynomialInFD<P> fftb;
    TwistIFFT<P>(fftb, b);
    MulInFD<P::n>(ffta, fftb);
    TwistFFT<P>(res, ffta);
}


template <class P>
inline void PolyMul(Polynomial<P> &res, const Polynomial<P> &a,
                    const Polynomial<P> &b)
{
    if constexpr (std::is_same_v<typename P::T, uint32_t>) {
        PolyMulFFT<P>(res, a, b);
    }
    else if constexpr (std::is_same_v<typename P::T, uint64_t>) {
        PolyMulNtt<P>(res, a, b);
    }
    else {
        PolyMulNaive<P>(res, a, b);
    }
}

template <class P>
inline void PolyMulRescaleUnsigned(Polynomial<P> &res,
                                   const UnsignedPolynomial<P> &a,
                                   const UnsignedPolynomial<P> &b)
{
    PolynomialInFD<P> ffta, fftb;
    TwistIFFT<P>(ffta, a);
    TwistIFFT<P>(fftb, b);
    MulInFD<P::n>(ffta, ffta, fftb);
    TwistFFTrescale<P>(res, ffta);
}

template <class P>
std::unique_ptr<std::array<PolynomialInFD<P>, 2 * P::n>> XaittGen()
{
    std::unique_ptr<std::array<PolynomialInFD<P>, 2 * P::n>> xaitt =
        std::make_unique<std::array<PolynomialInFD<P>, 2 * P::n>>();
    for (int i = 0; i < 2 * P::n; i++) {
        std::array<typename P::T, P::n> xai = {};
        xai[0] = -1;
        if (i < P::n)
            xai[i] += 1;
        else
            xai[i - P::n] -= 1;
        TwistIFFT<P>((*xaitt)[i], xai);
    }
    return xaitt;
}

template <class P>
std::unique_ptr<std::array<PolynomialNTT<P>, 2 * P::n>> XaittGenNTT()
{
    std::unique_ptr<std::array<PolynomialNTT<P>, 2 * P::n>> xaitt =
        std::make_unique<std::array<PolynomialNTT<P>, 2 * P::n>>();
    for (int i = 0; i < 2 * P::n; i++) {
        std::array<typename P::T, P::n> xai = {};
        xai[0] = -1;
        if (i < P::n)
            xai[i] += 1;
        else
            xai[i - P::n] -= 1;
        TwistINTT<P>((*xaitt[i]), xai);
    }
    return xaitt;
}

}  // namespace TFHEpp