// Independent enumeration by x and by the residue of h. No factorization of F.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;

std::ostream& write_wide(std::ostream& out, Wide value) {
    if (value >= 10) write_wide(out,value/10);
    return out << static_cast<int>(value%10);
}

Int integer_sqrt(Int n) {
    Int root = std::sqrt(static_cast<long double>(n));
    while ((root+1)*(root+1) <= n) ++root;
    while (root*root > n) --root;
    return root;
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int limit = std::atoi(argv[1]), cutoff = std::atoi(argv[2]);
    if (limit < 2 || limit > 4096 || cutoff < 0 || cutoff > 33554432) return 2;
    Int systems = 0, parity = 0;
    Wide largest = 0;
    for (Int m = 2; m <= limit; ++m) {
        Int square = m*m;
        for (Int u = 1; u < m; ++u) {
            Int k = m-u, km = k*m;
            Int last = (square+integer_sqrt(km*square-u))/u;
            Int first = std::max<Int>(cutoff+1, square/u+1);
            if (last < first) continue;
            std::vector<Int> divisors;
            for (Int d = 1; d*d <= km; ++d) {
                if (km%d) continue;
                divisors.push_back(d);
                if (d*d != km) divisors.push_back(km/d);
            }
            for (Int x = first; x <= last; ++x) {
                Int ell = u*x-square, top = square*(x-m);
                Int cap = std::min(km,top-ell*x);
                if (cap < 1) continue;
                auto check = [&](Int h) {
                    Int y = (top-h)/ell;
                    if (std::gcd(x,y) != 1 || h%std::gcd(m,x)) return;
                    Wide numerator = Wide(k)*x*y;
                    if (numerator%h) return;
                    Wide n = numerator/h, g = km/h, a = g*x, b = g*y;
                    ++systems;
                    largest = std::max(n,largest);
                    if ((a & (n-a)) && (b & (n-b))) return;
                    ++parity;
                    write_wide(std::cout,n) << ' ';
                    write_wide(std::cout,a) << ' ';
                    write_wide(std::cout,b) << '\n';
                };
                Int start = top%ell;
                if (!start) start = ell;
                if (start > cap) continue;
                if ((cap-start)/ell < static_cast<Int>(divisors.size())) {
                    for (Int h = start; h <= cap; h += ell) {
                        if (km%h == 0) check(h);
                    }
                } else {
                    for (Int h : divisors) {
                        if (h <= cap && (top-h)%ell == 0) check(h);
                    }
                }
            }
        }
    }
    std::cout << "Systems " << systems << " parity " << parity << " largest ";
    write_wide(std::cout,largest) << '\n';
}
