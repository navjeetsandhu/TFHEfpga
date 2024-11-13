#include <tfhe++.hpp>


template <class P>
void offsetGen() {
    constexpr typename P::T offset = TFHEpp::offsetgen<P>();
    std::cout << std::dec << "offset " << offset << " Ox"<< std::hex << offset << std::endl;
}

int main()
{
    offsetGen<TFHEpp::lvl1param>();
    offsetGen<TFHEpp::lvl2param>();
    offsetGen<TFHEpp::lvl3param>();

    return 0;
}
