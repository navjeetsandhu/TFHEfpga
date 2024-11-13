#include <tfhe++.hpp>
#include <chrono>
#include <iostream>
using namespace TFHEpp;
using namespace std;

uint8_t AndCheck(const uint8_t in0, const uint8_t in1) { return in0 * in1; }


void RunTest()
{
    cout << "------ Testing of FFT AND gate for lvl2param ------" << endl;
    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl02param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(*sk);
    vector<uint8_t> p = {1,1,0,0};
    vector<uint8_t> p2(4);
    vector<TLWE<lvl2param>> cres(3);
    vector<TLWE<lvl2param>> c(4);

    c = bootsSymEncrypt<lvl2param>(p, *sk);

    TFHEpp::HomAND<lvl20param, lvl02param, lvl2param::mu>(cres[0], c[0], c[1], ek);
    TFHEpp::HomAND<lvl20param, lvl02param, lvl2param::mu>(cres[1], c[0], c[2], ek);
    TFHEpp::HomAND<lvl20param, lvl02param, lvl2param::mu>(cres[2], c[3], c[2], ek);

    p2 = bootsSymDecrypt<lvl2param>(cres, *sk);

    c_assert(1 == p2[0]);
    c_assert(0 == p2[1]);
    c_assert(0 == p2[2]);
    cout << "------ pass ------" << endl;
}

int main()
{
    RunTest();
    return 0;
}