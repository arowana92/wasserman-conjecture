// Enumerate odd cofactors d in 2^h+1=dQ via the low binary digits of Q.
#include <cstdint>
#include <iostream>
#include <string>
using Word = std::uint64_t;
using Wide = __uint128_t;

struct OrderSearch {
    int bits;
    Word prefixes = 0, cofactors = 0, hits = 0, longest = 0;

    void complete(Word q, int previous, int carry) {
        ++prefixes;
        // Newton inversion doubles the number of correct low bits at every step.
        Word inverse = 1;
        for (int i = 0; i < 6; ++i) inverse *= 2 - q*inverse;
        Word d = inverse & ((Word(1) << bits)-1);
        if (d == 3 || d%24 != 3) return;
        ++cofactors;
        Word residual = (Wide(d)*q-1) >> bits;
        Word h = bits;
        // d*(known prefix) = 1 + 2^h * residual.
        while (residual != 1) {
            int bit = residual%2;
            int sum = bit + previous + carry;
            if (previous && sum%2) break;
            residual = (residual + bit*d)/2;
            previous = bit;
            carry = sum/2;
            ++h;
        }
        if (h > longest) longest = h;
        if (residual == 1 && h%2) {
            ++hits;
            std::cout << "hit " << d << ' ' << h << '\n';
        }
    }

    void visit(int index, Word q, int previous, int carry) {
        if (index == bits) {
            complete(q,previous,carry);
            return;
        }
        for (int bit = 0; bit <= 1; ++bit) {
            int sum = bit + previous + carry;
            if (!(previous && sum%2))
                visit(index+1,q+(Word(bit)<<index),bit,sum/2);
        }
    }
};

int main(int argc, char** argv) {
    int bits = argc == 2 ? std::stoi(argv[1]) : 42;
    if (argc > 2 || bits < 3 || bits > 55) return 2;
    OrderSearch search{bits};
    search.visit(1,1,1,0);
    std::cout << "prefixes, cofactors, hits, last position\n";
    std::cout << search.prefixes << ' ' << search.cofactors << ' '
              << search.hits << ' ' << search.longest << '\n';
}
