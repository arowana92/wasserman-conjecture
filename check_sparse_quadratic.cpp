// At full power P=p only the quadratic branch lies above 5p.
#include "range_factors.hpp"
#include <cstdint>
#include <cstdlib>
#include <iostream>

using U64 = uint64_t;
using Wide = __uint128_t;

bool avoids(U64 singleton, U64 repeated, int base) {
    while (singleton || repeated) {
        if (singleton % base + 2*(repeated % base) >= U64(base)) return false;
        singleton /= base;
        repeated /= base;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int limit = std::atoi(argv[1]);
    if (limit < 5 || limit > 200000000) return 2;
    FactorSieve sieve(limit + 1);
    U64 primes = 0, nontrivial = 0, feasible = 0;
    for (int p = 5; p <= limit; p += 2) {
        if (!sieve.prime(p)) continue;
        ++primes;
        U64 singleton = U64(p)*p, repeated = singleton+p;
        U64 n = singleton+2*repeated, S = 1;
        for (auto [q, w] : sieve.factor(p+1)) {
            // At 2 the two identical positive parts always cause a carry.
            if (q != 2 && avoids(singleton, repeated, q)) S *= w;
        }
        if (S > 1) ++nontrivial;
        U64 T = (n-1)/S;
        if (Wide(27)*T*T > Wide(p)*p*p) continue;
        std::cout << "CRT " << p << " " << p << " " << n << " " << S << "\n";
        ++feasible;
    }
    std::cout << "Primes " << primes << ", nontrivial predecessor products "
              << nontrivial << ", passing the companion size bound " << feasible << "\n";
}
