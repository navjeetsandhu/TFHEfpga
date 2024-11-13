#include "c_assert.hpp"
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    cout << "test p=1" << endl;

    cout << "lvl1" << endl;
    lweKey key;

    array<bool, lvl1param::n> p;
    for (bool &i : p) i = (binary(engine) > 0);
    Polynomial<lvl1param> pmu;
    for (int i = 0; i < lvl1param::n; i++)
        pmu[i] = p[i] ? lvl1param::mu : -lvl1param::mu;
    TRLWE<lvl1param> c = trlweSymEncrypt<lvl1param>(pmu, key.lvl1);

    const Polynomial<TFHEpp::lvl1param> plainpoly = {
        static_cast<typename lvl1param::T>(1)};

    TRGSWFFT<lvl1param> trgswfft =
        trgswfftSymEncrypt<lvl1param>(plainpoly, key.lvl1);
    trgswfftExternalProduct<lvl1param>(c, c, trgswfft);
    array<bool, lvl1param::n> p2 = trlweSymDecrypt<lvl1param>(c, key.lvl1);
    for (int i = 0; i < lvl1param::n; i++) {
        c_assert(p[i] == p2[i]);
    }

    cout << "Passed" << endl;
}