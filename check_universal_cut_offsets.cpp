// Enumerate the possible weights independently using V = U^{-1} modulo S.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

using Int = std::int64_t;
using Factors = std::vector<std::pair<Int,Int>>;

Factors factor(Int n) {
    Factors result;
    for (Int p = 2; p <= n/p; ++p) {
        if (n % p) continue;
        Int power = 1;
        do { n /= p; power *= p; } while (n % p == 0);
        result.emplace_back(p, power);
    }
    if (n > 1) result.emplace_back(n, n);
    return result;
}

std::vector<Int> divisors(Int n) {
    std::vector<Int> result{1};
    for (auto [p, power] : factor(n)) {
        auto count = result.size();
        for (Int w = p; w <= power; w *= p)
            for (size_t i = 0; i < count; ++i)
                result.push_back(result[i]*w);
    }
    return result;
}

Int inverse(Int a, Int modulus) {
    Int x = 1, y = 0, b = modulus;
    while (b) {
        Int q = a/b;
        a -= q*b;
        x -= q*y;
        std::swap(a,b);
        std::swap(x,y);
    }
    if (a != 1) throw std::logic_error("inverse requires coprime factors");
    return (x % modulus + modulus) % modulus;
}

// Compare factorial floors at every layer, allowing a zero third argument for a cut.
bool unit(Int n, Int a, Int b, Int c, Int p) {
    while (n >= p) {
        n /= p; a /= p; b /= p; c /= p;
        if (n != a+b+c) return false;
    }
    return true;
}

bool avoids_factors(Int n, Int a, Int b, Int c, const Factors& factors) {
    for (auto [p,w] : factors)
        if (!unit(n,a,b,c,p)) return false;
    return true;
}

using Hit = std::array<Int,5>;

void companions(Int n, Int a, Int u, Int v, Int t, std::set<Hit>& hits) {
    Factors shifted = factor(t);
    Int largest = 1;
    for (auto [p,w] : shifted) largest = std::max(largest,w);
    if (2*largest >= u) return;
    Factors weight = factor(v);
    // Some pair of normalized parts are positive multiples of 'largest'.
    // When t=1, this just enumerates the complete set of positive triples.
    for (Int i = largest; 2*i < u; i += largest) {
        for (Int j = i; i+j < u; j += largest) {
            Int x=v*i, y=v*j, z=n-x-y;
            if (!avoids_factors(n,x,y,z,shifted) || !avoids_factors(n,x,y,z,weight))
                continue;
            std::array<Int,3> row{x,y,z};
            std::sort(row.begin(),row.end());
            hits.insert({n,a,row[0],row[1],row[2]});
        }
    }
}

void search(Int denominator, Int radius) {
    std::set<Hit> hits;
    for (Int d = -radius; d <= radius; ++d) {
        if (!d || d == -1 || d == denominator-1 || (denominator == 2 && d > 0)) continue;
        Int lcm = std::lcm(std::abs(d+1), std::abs(d+1-denominator));
        auto possible_s = divisors(lcm);
        for (Int u : divisors(std::abs(d))) {
            if (u < 3) continue;
            auto cut_weight = factor(u);
            for (Int s : possible_s) {
                if (std::gcd(u,s) != 1) continue;
                auto cut_shift = factor(s);
                Int first = s == 1 ? 1 : inverse(u,s);
                for (Int v = first; ; v += s) {
                    Int n = u*v, t = (n-1)/s;
                    if (27*t*t > u*u*u) break;
                    if (v < 2 || (n+d) % denominator || std::gcd(u,v)>1 || std::gcd(s,t)>1)
                        continue;
                    Int a = (n+d)/denominator;
                    if (a <= 0 || a >= n || a % u) continue;
                    if (avoids_factors(n,a,n-a,0,cut_weight) &&
                        avoids_factors(n,a,n-a,0,cut_shift))
                        companions(n,a,u,v,t,hits);
                }
            }
        }
    }
    std::cout << "denominator " << denominator << ", radius " << radius
              << ", distinct cut-row witnesses " << hits.size() << '\n';
    for (auto hit : hits) {
        for (auto x : hit) std::cout << x << ' ';
        std::cout << '\n';
    }
}

int main() {
    for (auto [t,h] : std::array<std::pair<Int,Int>,5>{{
             {2,4616},{3,4048},{4,113},{5,2912},{6,2344}}})
        search(t,h);
}
