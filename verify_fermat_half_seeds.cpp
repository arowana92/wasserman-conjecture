// Enumerate z by its own low bits and lift (P*z)^-1 by Newton's identity.
#include <array>
#include <cstdint>
#include <iostream>
#include <string>

using Wide = unsigned __int128;
using Word = std::uint64_t;
constexpr Word prime = 65537;
constexpr Word x = 49153;
constexpr Word y = 16384;
constexpr Word odd_modulus = 13 * prime;

std::string decimal(Wide value) {
    if (value < 10) return std::to_string(static_cast<unsigned>(value));
    return decimal(value / 10) + char('0' + value % 10);
}

struct Search {
    unsigned b, width;
    Wide modulus, limit;
    Word inverse_block = 1;
    std::array<Word, 5> residues{};

    explicit Search(unsigned exponent)
        : b(exponent), width((b + 15) / 2),
          modulus(Wide(1) << (b + 14)), limit((modulus - 1) / prime) {
        if (width < 15) width = 15;
        // Divide by 2 modulo 13*P, once for each digit of the prefix.
        for (unsigned i = 0; i < width; ++i) {
            if (inverse_block % 2) inverse_block += odd_modulus;
            inverse_block /= 2;
        }
        Word r_unit = static_cast<Word>((Wide(1) << (b-2)) % prime);
        for (Word digit = 0; digit < residues.size(); ++digit) {
            Word r = digit * r_unit % prime;
            while (r % 13) r += prime;
            residues[digit] = r;
        }
    }

    void visit(unsigned length, Word z, Word t) const {
        if (length < width) {
            Word place = Word(1) << length;
            Word mask = 2*place - 1;
            for (Word next_z : {z, z + place}) {
                Word next_t = t;
                if ((prime*next_z*t - 1) & place) next_t += place;
                Word sum = (x*next_t & mask) + (y*next_t & mask);
                if (sum <= mask) visit(length+1, next_z, next_t);
            }
            return;
        }
        Wide block = Wide(1) << width;
        Wide mask = modulus-1;
        for (Word residue : residues) {
            Word offset = (residue + odd_modulus - z % odd_modulus) % odd_modulus;
            Word index = offset * inverse_block % odd_modulus;
            for (Wide part = z + block*index; part <= limit; part += block*odd_modulus) {
                Wide scale = Wide(t)*(2 - Wide(prime)*part*t) & mask;
                // Check every low digit. An AND detects a carry even if it has
                // disappeared again by the top of the modulus.
                if ((x*scale & y*scale & mask) == 0)
                    std::cout << b << ' ' << decimal(part) << '\n';
            }
        }
    }
};

int main() {
    for (unsigned b = 4; b <= 68; b += 4)
        Search(b).visit(15, 16385, 16385);
}
