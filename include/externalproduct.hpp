#pragma once

#include <array>
#include <cstdint>
#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "decomposition.hpp"


namespace TFHEpp {

template <class P>
void trgswfftExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWFFT<P> &trgswfft)
{
    alignas(64) DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, trlwe[0]);
    alignas(64) PolynomialInFD<P> decpolyfft;
    TwistIFFT<P>(decpolyfft, decpoly[0]);
    alignas(64) TRLWEInFD<P> restrlwefft;
    for (int m = 0; m < P::k + 1; m++)
        MulInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[0][m]);
    for (int i = 1; i < P::l; i++) {
        TwistIFFT<P>(decpolyfft, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i][m]);
    }
    for (int k = 1; k < P::k + 1; k++) {
        Decomposition<P>(decpoly, trlwe[k]);
        for (int i = 0; i < P::l; i++) {
            TwistIFFT<P>(decpolyfft, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
                FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                              trgswfft[i + k * P::l][m]);
        }
    }
    for (int k = 0; k < P::k + 1; k++) TwistFFT<P>(res[k], restrlwefft[k]);
}

template <class P>
void trgswnttExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWNTT<P> &trgswntt)
{
    DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, trlwe[0]);
    PolynomialNTT<P> decpolyntt;
    TwistINTT<P>(decpolyntt, decpoly[0]);
    TRLWENTT<P> restrlwentt;
    for (int m = 0; m < P::k + 1; m++)
        for (int i = 0; i < P::n; i++)
            restrlwentt[m][i] = decpolyntt[i] * trgswntt[0][m][i];
    for (int i = 1; i < P::l; i++) {
        TwistINTT<P>(decpolyntt, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            for (int j = 0; j < P::n; j++)
                restrlwentt[m][j] += decpolyntt[j] * trgswntt[i][m][j];
    }
    for (int k = 1; k < P::k + 1; k++) {
        Decomposition<P>(decpoly, trlwe[k]);
        for (int i = 0; i < P::l; i++) {
            TwistINTT<P>(decpolyntt, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] +=
                        decpolyntt[j] * trgswntt[i + k * P::l][m][j];
        }
    }
    for (int k = 0; k < P::k + 1; k++) TwistNTT<P>(res[k], restrlwentt[k]);
}

template <class P>
void trgswrainttExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                                const TRGSWRAINTT<P> &trgswntt)
{
    DecomposedPolynomialRAINTT<P> decpolyntt;
    DecompositionRAINTT<P>(decpolyntt, trlwe[0]);
    TRLWERAINTT<P> restrlwentt;
    for (int m = 0; m < P::k + 1; m++)
        for (int i = 0; i < P::n; i++)
            restrlwentt[m][i] =
                raintt::MulSREDC(decpolyntt[0][i], trgswntt[0][m][i]);
    for (int i = 1; i < P::l; i++) {
        for (int m = 0; m < P::k + 1; m++)
            for (int j = 0; j < P::n; j++)
                restrlwentt[m][j] = raintt::AddMod(
                    restrlwentt[m][j],
                    raintt::MulSREDC(decpolyntt[i][j], trgswntt[i][m][j]));
    }
    for (int k = 1; k < P::k + 1; k++) {
        DecompositionRAINTT<P>(decpolyntt, trlwe[k]);
        for (int i = 0; i < P::l; i++) {
            for (int m = 0; m < P::k + 1; m++)
                for (int j = 0; j < P::n; j++)
                    restrlwentt[m][j] = raintt::AddMod(
                        restrlwentt[m][j],
                        raintt::MulSREDC(decpolyntt[i][j],
                                         trgswntt[i + k * P::l][m][j]));
        }
    }
    // if constexpr(hasq<P>) for (int k = 0; k < P::k + 1; k++)
    // raintt::TwistNTT<typename P::T,P::nbit, P::q == (1ULL<<P::qbit)>(res[k],
    // restrlwentt[k],(*raintttable)[0],(*raintttwist)[0]);
    for (int k = 0; k < P::k + 1; k++)
        raintt::TwistNTT<typename P::T, P::nbit, true>(
            res[k], restrlwentt[k], (*raintttable)[0], (*raintttwist)[0]);
}

}  // namespace TFHEpp