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

    cout << "test p=-1" << endl;
    cout << "lvl2" << endl;

    lweKey key;

    array<bool, lvl2param::n> p;
    for (bool &i : p) i = binary(engine) > 0;
    Polynomial<lvl2param> pmu;
    for (int i = 0; i < lvl2param::n; i++)
        pmu[i] = p[i] ? lvl2param::mu : -lvl2param::mu;
    TRLWE<lvl2param> c = trlweSymEncrypt<lvl2param>(pmu, key.lvl2);

    const Polynomial<TFHEpp::lvl2param> plainpoly = {
        static_cast<typename lvl2param::T>(-1)};

    TRGSWFFT<lvl2param> trgswfft =
        trgswfftSymEncrypt<lvl2param>(plainpoly, key.lvl2);
    trgswfftExternalProduct<lvl2param>(c, c, trgswfft);
    array<bool, lvl2param::n> p2 = trlweSymDecrypt<lvl2param>(c, key.lvl2);
    for (int i = 0; i < lvl2param::n; i++) {
        c_assert(p[i] == !p2[i]);
    }

    cout << "Passed" << endl;
}