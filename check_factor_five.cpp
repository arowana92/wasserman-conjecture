// Generate every binary and CRT candidate for the predecessor power five.
// Each line of stdout is the odd factor, the binary exponent and the CRT multiplier.
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

using Word = std::uint64_t;

struct Search {
    int bits, max_k;
    Word shapes = 0, compatible = 0, candidates = 0;

    void inspect(Word q) {
        ++shapes;
        if (q <= 5) return;
        int first = 3;
        while (first <= 6 && ((q % 15) << first) % 15 != 7) ++first;
        if (first == 7 || first > max_k) return;
        ++compatible;
        // 5*2^first is at most 320. Solve for its inverse using the small modulus.
        Word m = Word(5) << first;
        Word digit = 0;
        while ((digit*(q % m) + 1) % m) ++digit;
        Word inverse = (digit*q + 1)/m;
        Word correction = 1;
        while (((q % 16)*correction + 1) % 16) ++correction;
        constexpr Word pow2mod9[] = {1, 2, 4, 8, 7, 5};
        for (int k = first; k <= max_k; k += 4) {
            Word l = inverse;
            inverse = (inverse + q*((inverse*correction) % 16))/16;
            Word z = 5*l;
            if (l % 3 || (z & ~(2*q))) continue;
            if ((5*(q % 9)*pow2mod9[k % 6]+1) % 9 == 0) continue;
            ++candidates;
            std::cout << q << ' ' << k << ' ' << l << '\n';
        }
    }

    // Add q and 2q from low to high. The sum 3q must avoid the bits of 2q.
    void visit(int place, Word q, int previous, int carry) {
        if (place == bits) {
            if (!(previous && (previous+carry) % 2)) inspect(q);
            return;
        }
        for (int bit = 0; bit <= 1; ++bit) {
            int total = bit + previous + carry;
            if (previous && total % 2) continue;
            visit(place+1, q+(Word(bit) << place), bit, total/2);
        }
    }
};

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int bits = std::stoi(argv[1]), max_k = std::stoi(argv[2]);
    // Products digit*q in the first inverse and 15*q in its update fit in 64 bits.
    if (bits < 1 || bits > 54 || max_k < 3 || max_k > 1000000) return 2;
    Search search{bits, max_k};
    search.visit(1, 1, 1, 0);
    std::cerr << "Binary shapes " << search.shapes
              << ", compatible residue classes " << search.compatible
              << ", candidate triples " << search.candidates << '\n';
}
