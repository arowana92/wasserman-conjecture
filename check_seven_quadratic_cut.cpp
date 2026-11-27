// Independent search: choose the smallest CRT group before choosing P.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>

using Int = std::int64_t;
struct Totals { Int cuts = 0, parity = 0, rows = 0; };

Int residue(Int x, Int n) {
    Int r = x%n;
    return r < 0 ? r+n : r;
}

Int inv(Int x, Int n) {
    if (n == 1) return 0;
    if (x == 1) return 1;
    // Invert n mod x recursively, then solve x*y = 1 mod n.
    return residue((1-n*inv(n%x,x))/x,n);
}

Totals scan(int limit, int t) {
    const int u = 7-t;
    Totals totals;
    for (Int d = 3; d*d*d <= 10LL*limit+7; d += 2) {
        if (std::gcd(d,10LL) != 1) continue;
        const Int step = d*d;
        const Int ten_inverse = inv(10%step,step);
        // All normalized numerators at a nonempty smallest group are positive.
        for (Int n = 1; t*n < d; ++n) {
            if (std::gcd(n,d) != 1) continue;
            Int m0 = residue(-u*inv(t*n%d,d),d);
            Int first = residue((d*m0-7)*ten_inverse,step);
            if (first == 0) first = step;
            for (Int P = first; P <= limit; P += step) {
                if (t == 5 && P%2 == 0 && (P & (P-1))) continue;
                Int H = 10*P+7;
                Int m = H/d;
                if (m < d*d) continue;
                Int W = 1+u*P, T = 1+t*P, N = W*T;
                Int f = m*n, Y = f*T;
                if (f >= W) continue;
                ++totals.cuts;
                bool need_odd = t == 2 || P%2 == 1;
                if (need_odd && (Y & (N-Y))) continue;
                ++totals.parity;
                // Split the remaining integer by ordinary trial division;
                // neither its factorization nor a prime sieve is needed.
                for (Int a = d; a*a <= m; ++a) {
                    if (m%a) continue;
                    Int b = m/a;
                    if (std::gcd(a,b) != 1) continue;
                    Int other = d*b;
                    Int f2 = other*(residue(W,a)*inv(other%a,a)%a);
                    Int f3 = W-f-f2;
                    if (f3 <= 0) continue;
                    std::array<Int,3> C = {Y,T*f2,T*f3};
                    if (need_odd && (C[1] & C[2])) continue;
                    ++totals.rows;
                    std::cout << "Unexcluded P=" << P << " t=" << t << " companion "
                              << C[0] << " " << C[1] << " " << C[2] << "\n";
                }
            }
        }
    }
    return totals;
}

int main(int argc, char** argv) {
    int limit = argc == 2 ? std::stoi(argv[1]) : 10000000;
    if (limit < 1 || limit > 100000000) throw std::invalid_argument("bound must be in [1,100000000]");
    Int survivors = 0;
    for (int t : {2,5}) {
        auto c = scan(limit,t);
        std::cout << "P through " << limit << ", t=" << t
                  << ", cuts: " << c.cuts << ", after parity: " << c.parity
                  << ", complete rows after parity: " << c.rows << "\n";
        survivors += c.rows;
    }
    return survivors != 0;
}
