// List necessary row and exponent data at d*2^a+1 for odd d through 65535.
// The independent Python enumeration starts with x+y instead of d.
#include "range_factors.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

using Integer = std::int64_t;

Integer residue(Integer n, Integer m) {
    return (n % m + m) % m;
}

Integer inverse(Integer n, Integer m) {
    Integer a = n, b = m, x = 1, y = 0;
    while (b) {
        Integer quotient = a / b;
        std::tie(a, b) = std::make_tuple(b, a - quotient*b);
        std::tie(x, y) = std::make_tuple(y, x - quotient*y);
    }
    return residue(x, m);
}

struct Progression {
    Integer first = 0, period = 1;
};

class BinaryPowers {
    std::map<int, std::vector<int>> tables;
public:
    // Intersect the present exponent progression with 2^a = target (mod m).
    bool restrict(Progression& exponents, int m, int target) {
        auto& log = tables[m];
        if (log.empty()) {
            log.assign(m+1, -1);
            int r = 1, order = 0;
            do {
                log[r] = order++;
                r = 2*r % m;
            } while (r != 1);
            log[m] = order;
        }
        if (log[target] < 0) return false;
        Integer diff = log[target] - exponents.first;
        Integer g = std::gcd<Integer>(exponents.period, log[m]);
        if (diff % g) return false;
        Integer k = log[m]/g;
        Integer jump = residue(diff/g, k) * inverse(exponents.period/g % k, k) % k;
        exponents.first += jump*exponents.period;
        exponents.period *= k;
        return true;
    }
};

bool companion_fits(const Factors& factors, int u, int index, int modulus,
                    std::array<int, 3> coordinates) {
    int minimum = 0;
    for (int r : coordinates) minimum += r ? r : modulus;
    if (minimum > u) return false;
    if (index == factors.size) return true;
    int w = factors.atoms[index].power;
    Integer inv = inverse(modulus % w, w);
    for (int special = 0; special < 3; ++special) {
        auto next = coordinates;
        for (int j = 0; j < 3; ++j) {
            int target = j == special ? u % w : 0;
            next[j] += modulus * (residue(target-next[j], w)*inv % w);
        }
        if (companion_fits(factors, u, index+1, modulus*w, next)) return true;
    }
    return false;
}

void assign(const FactorSieve& sieve, BinaryPowers& powers, const Factors& atoms,
            int d, int x, int y, int u, int index, int t, Progression exponents) {
    if (index == atoms.size) {
        if (companion_fits(sieve.factor(t), u, 0, 1, {0, 0, 0})) {
            std::cout << d << ' ' << x << ' ' << y << ' ' << d-x-y << ' '
                      << u << ' ' << t << ' ' << exponents.first << ' '
                      << exponents.period << '\n';
        }
        return;
    }
    auto [p, w] = atoms.atoms[index];
    int next_t = t*w;
    if (2*w+1 <= u && Integer(27)*next_t*next_t <= Integer(u)*u*u)
        assign(sieve, powers, atoms, d, x, y, u, index+1, next_t, exponents);

    if (x % w && y % w) return;
    int other = x % w ? x : y;
    if (other % w && (other % p == 0 || !powers.restrict(exponents, w, inverse(other % w, w))))
        return;
    assign(sieve, powers, atoms, d, x, y, u, index+1, t, exponents);
}

int main(int argc, char** argv) {
    int limit = argc == 2 ? std::stoi(argv[1]) : 65535;
    if (argc > 2 || limit < 1 || limit > 65535) return 2;
    FactorSieve sieve(limit);
    BinaryPowers powers;
    std::vector<std::vector<int>> divisors(limit+1);
    for (int u = 3; u <= limit; u += 2)
        for (int k = u; k <= limit; k += u) divisors[k].push_back(u);

    for (int d = 3; d <= limit; d += 2) {
        Factors atoms = sieve.factor(d);
        if (atoms.size < 2) continue;  // N-1 then has at most two prime divisors.
        for (int x = (d-1)&d; x; x = (x-1)&d) {
            int rest = d^x;
            for (int y = (rest-1)&rest; y > x; y = (y-1)&rest) {
                int g = std::gcd(x, y);
                if (divisors[g].empty()) continue;
                // Each predecessor power assigned to the first row divides xy.
                Integer required = d/std::gcd<Integer>(d, Integer(x)*y);
                Integer minimum = 27*required*required;
                if (minimum > Integer(g)*g*g) continue;
                for (int u : divisors[g]) {
                    if (std::gcd(d, u) != 1 || minimum > Integer(u)*u*u) continue;
                    Progression exponents;
                    if (powers.restrict(exponents, u, u-inverse(d % u, u)))
                        assign(sieve, powers, atoms, d, x, y, u, 0, 1, exponents);
                }
            }
        }
    }
}
