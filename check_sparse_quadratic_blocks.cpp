// Check the binary cofactor theorem using actual integer weights at each scale.
#include "range_factors.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using U64 = std::uint64_t;
constexpr int limit = 1000000;
constexpr int known_scale = 200000000;

int inverse(U64 value, int modulus) {
    std::int64_t a = value % modulus, b = modulus, x = 1, y = 0;
    while (b) {
        auto quotient = a/b;
        auto rest = a-quotient*b, next = x-quotient*y;
        a = b;
        b = rest;
        x = y;
        y = next;
    }
    assert(a == 1);
    return (x % modulus + modulus) % modulus;
}

int main() {
    FactorSieve sieve(limit);
    U64 singles = 0, pairs = 0, finite = 0;
    for (int d = 1; d <= limit; d += 2) {
        if (d % 3 == 0) continue;
        U64 q = 2;
        while (q <= 2U*d || (d*q+1) % 3 != 0) q *= 2;
        U64 v = d*q+3, total = d*((d*q+4)/3);
        // Pick just the low ten bits of v*i as a submask, then visit every
        // positive i in the corresponding residue class.
        unsigned mask = total % 1024;
        unsigned inverse_v = inverse(v, 1024);
        std::vector<unsigned> entries;
        for (unsigned low = mask;; low = (low-1) & mask) {
            unsigned first = (low*inverse_v) % 1024;
            for (unsigned i = first ? first : 1024; 3*i < unsigned(d); i += 1024) {
                if ((v*i & ~total) == 0) entries.push_back(i);
            }
            if (!low) break;
        }
        singles += entries.size();
        for (unsigned a = 0; a < entries.size(); ++a) {
            for (unsigned b = a+1; b < entries.size(); ++b) {
                U64 i = entries[a], j = entries[b];
                if (3*(i+j) >= unsigned(d) || (v*i & v*j)) continue;
                ++pairs;
                assert(i*j % d != 0);
            }
        }

        // Choose the divisor of d on which 3*q*i must be one. On its
        // complementary divisor i vanishes, fixing i uniquely below d.
        std::vector<int> divisors{1};
        for (auto atom : sieve.factor(d)) {
            auto size = divisors.size();
            for (unsigned j = 0; j < size; ++j) divisors.push_back(divisors[j]*atom.power);
        }
        for (q /= 4; q >= 2 && d*q+1 > 3U*known_scale; q /= 4) {
            assert(q <= 2U*d && (d*q+1) % 3 == 0);
            ++finite;
            v = d*q+3;
            total = d*((d*q+4)/3);
            for (int b : divisors) {
                if (b == 1) continue;
                int a = d/b;
                U64 i = U64(a)*inverse(3*q*a, b);
                if (3*i < unsigned(d)) {
                    U64 part = v*i;
                    assert(__builtin_popcountll(part) + __builtin_popcountll(total-part)
                           > __builtin_popcountll(total));
                }
            }
        }
    }
    assert(singles == 15433 && pairs == 3 && finite == 1471770);
    std::cout << singles << " entries, " << pairs << " pairs, " << finite
              << " finite scales, all excluded.\n";
}
