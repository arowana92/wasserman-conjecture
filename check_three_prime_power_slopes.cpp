// c++ -std=c++17 -O3 check_three_prime_power_slopes.cpp -o check_slopes
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;

Int inverse(Int a, Int m) {
    Int b = m, u = 1, v = 0;
    while (b) {
        Int k = a / b;
        Int rem = a % b, next = u - k * v;
        a = b; b = rem; u = v; v = next;
    }
    assert(a == 1);
    return (u % m + m) % m;
}

bool subword(Int part, Int total, Int base) {
    if (base == 2) return (part & ~total) == 0;
    for (; part; part /= base, total /= base)
        if (part % base > total % base) return false;
    return true;
}

// Generate multiples of divisor that can occupy some digits of total.
std::vector<Int> multipliers(Int total, Int divisor, Int base, Int own_prime) {
    std::vector<std::pair<Int, Int>> digits;
    for (Int n = total, place = 1; n; n /= base, place *= base)
        if (n % base) digits.emplace_back(n % base, place);
    std::sort(digits.rbegin(), digits.rend());
    std::vector<Int> left{0}, right{0};
    for (auto [digit, place] : digits) {
        auto &side = left.size() < right.size() ? left : right;
        auto old_size = side.size();
        for (Int d = 1; d <= digit; ++d)
            for (std::size_t i = 0; i < old_size; ++i)
                side.push_back(side[i] + d * place);
    }
    std::unordered_multimap<Int, Int> lookup;
    for (Int value : left) lookup.emplace(value % divisor, value);
    std::vector<Int> out;
    for (Int high : right) {
        auto range = lookup.equal_range((divisor - high % divisor) % divisor);
        for (auto it = range.first; it != range.second; ++it) {
            Int slope = high + it->second, x = slope / divisor;
            if (slope && slope < total && x % own_prime) out.push_back(x);
        }
    }
    return out;
}

Int mod_power(Int x, Int e, Int m) {
    Int out = 1;
    for (; e; e /= 2) {
        if (e % 2) out = Wide(out) * x % m;
        x = Wide(x) * x % m;
    }
    return out;
}

// All lifts of one exponent class must pass every floor layer fixed by it.
bool fixed_layers(Int p, Int r, Int period, Int exponent, Int first, Int second,
                  Int total, Int required_modulus, Int required_slope, Int residue) {
    for (Int m = p; mod_power(r, period, m) == 1; m *= p) {
        Int z = mod_power(r, exponent, m), d = std::min(m, required_modulus);
        if ((Wide(required_slope) * z + d - residue) % d != 0) return false;
        Int x = Wide(first) * z % m, y = Wide(second) * z % m;
        Int w = (Wide(total - first - second) * z + 1) % m;
        if (x + y + w >= m) return false;
        assert(m <= std::numeric_limits<Int>::max() / p);
    }
    return true;
}

bool possible_scale(Int p, Int q, Int r, Int P, Int Q, Int x, Int y, Int g,
                    Int period) {
    for (Int ell = 2; g > 1; ++ell) {
        if (ell * ell > g) ell = g;
        if (g % ell) continue;
        Int required = ell == 2 ? 4 : ell;
        bool can_assign = g % required == 0;
        do { g /= ell; } while (g % ell == 0);
        if (!can_assign) continue;
        Int u = x * Q, v = y * P, total = P * Q;
        for (Int e = 0; e < period; ++e) {
            if (fixed_layers(p, r, period, e, u, v, total, P, u, 1) &&
                fixed_layers(q, r, period, e, u, v, total, Q, v, 1) &&
                fixed_layers(ell, r, period, e, u, v, total, required, total, -1))
                return true;
        }
    }
    return false;
}

struct Counts { Int cells{}, slopes{}, first{}, second{}, remaining{}; };

void cell(Int p, Int q, Int r, Int P, Int Q, Int period, Counts &counts) {
    ++counts.cells;
    Int total = P * Q;
    auto X = multipliers(total, Q, r, p);
    if (X.empty()) return;
    auto Y = multipliers(total, P, r, q);
    if (Y.empty()) return;
    for (Int x : X) {
        Int ip = inverse(x, P), rp = ip * inverse(Q % P, P) % P;
        for (Int y : Y) {
            if (!subword(y * P, total - x * Q, r)) continue;
            Int g = std::gcd(x, y);
            while (g % r == 0) g /= r;
            if (g < 3) continue;
            ++counts.slopes;
            if (!subword(y * rp % P, ip, p)) continue;
            ++counts.first;
            Int iq = inverse(y, Q), rq = iq * inverse(P % Q, Q) % Q;
            if (!subword(x * rq % Q, iq, q)) continue;
            ++counts.second;
            if (possible_scale(p, q, r, P, Q, x, y, g, period)) {
                ++counts.remaining;
                std::cout << "survivor " << P << ' ' << Q << ' ' << r << ' '
                          << x << ' ' << y << '\n';
            }
        }
    }
}

int main() {
    const std::array<std::array<Int, 5>, 6> sets{{
        {2,3,5,10000000,12}, {2,3,7,10000000,12}, {2,5,7,10000000,12},
        {3,5,7,10000000,12}, {2,3,11,1000000,60}, {2,3,13,1000000,60}
    }};
    Int cells = 0, survivors = 0;
    for (auto set : sets) {
        Counts counts;
        for (int k = 0; k < 3; ++k) {
            Int p = std::min(set[(k+1)%3], set[(k+2)%3]);
            Int q = std::max(set[(k+1)%3], set[(k+2)%3]);
            for (Int P = p; P <= set[3]; P *= p)
                for (Int Q = q; Q <= set[3]; Q *= q)
                    cell(p, q, set[k], P, Q, set[4], counts);
        }
        std::cout << set[0] << ' ' << set[1] << ' ' << set[2] << " counts "
                  << counts.cells << ' ' << counts.slopes << ' ' << counts.first
                  << ' ' << counts.second << ' ' << counts.remaining << std::endl;
        cells += counts.cells;
        survivors += counts.remaining;
    }
    return cells != 2902 || survivors != 0;
}
