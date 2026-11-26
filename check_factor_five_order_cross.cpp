// Lift the cofactor bits directly, retaining every possible inverse prefix.
#include <cstdint>
#include <iostream>
#include <string>
using Word=std::uint64_t;
using Wide=__uint128_t;

struct Lift {
    int bits;
    Word prefixes=0, cofactors=0, long_words=0;
    void visit(int length, Word d, Word q, Word quotient) {
        if (length==bits) {
            ++prefixes;
            if (d==3 || d%24!=3) return;
            ++cofactors;
            Wide inverse=q;
            for (int i=length; i<128; i*=2) inverse*=2-Wide(d)*inverse;
            if (!((3*inverse)&(2*inverse))) {
                ++long_words;
                std::cout << "long " << d << '\n';
            }
            return;
        }
        Word place=Word(1)<<length;
        for (Word x=0; x<2; ++x) {
            Word y=(quotient+x)&1;
            Word next=q+y*place;
            if ((3*next)&(2*next)&place) continue;
            // (d+x*place)(q+y*place)-1, divided by twice the place.
            Word new_quotient=(quotient+x*q+y*d+x*y*place)/2;
            visit(length+1,d+x*place,next,new_quotient);
        }
    }
};

int main(int argc,char** argv) {
    int bits=argc==2 ? std::stoi(argv[1]):42;
    if (argc>2 || bits<3 || bits>55) return 2;
    Lift search{bits};
    search.visit(1,1,1,0);
    std::cout << "prefixes, cofactors, long words\n"
              << search.prefixes << ' ' << search.cofactors << ' ' << search.long_words << '\n';
}
