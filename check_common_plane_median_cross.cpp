// Independent search using the full prime powers of the largest coefficient.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

using Int = std::int64_t;
using Row = std::array<Int, 3>;
constexpr int bound = 211;
constexpr Int cutoff = 1000000000;

Int floor_sqrt(Int n) {
    Int r = static_cast<Int>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

// For each full prime power, an idempotent is either 0 or 1 modulo that power.
std::vector<int> idempotents(int t) {
    std::vector<int> residues{0};
    int remaining = t, modulus = 1;
    for (int p = 2; remaining > 1; ++p) {
        if (p > remaining / p) p = remaining;
        if (remaining % p != 0) continue;
        int q = 1;
        do {
            remaining /= p;
            q *= p;
        } while (remaining % p == 0);
        int basis = modulus;
        while (basis % q != 1) basis += modulus;
        std::vector<int> next;
        for (int old : residues) {
            Int zero = (Int(old) * (1 - basis)) % (modulus * q);
            if (zero < 0) zero += modulus * q;
            next.push_back(static_cast<int>(zero));
            next.push_back(static_cast<int>((zero + basis) % (modulus * q)));
        }
        residues.swap(next);
        modulus *= q;
    }
    return residues;
}

bool make_row(Int n, Int cut, int r, int s, int t, Row& row) {
    Int numerator = (t + r) * cut - r * n;
    if (numerator % (s - r)) return false;
    Int y = numerator / (s - r);
    row = {n - cut - y, y, cut};
    return row[0] > 0 && row[1] > 0;
}

Int factorial_twos(Int n) {
    Int value = 0;
    while (n) {
        n /= 2;
        value += n;
    }
    return value;
}

Int valuation(Int n, const Row& row) {
    Int v = factorial_twos(n);
    for (Int part : row) v -= factorial_twos(part);
    return v;
}

// Scan every possible smaller cut at each of the two reported weaker survivors.
void check_examples() {
    for (auto [r, s, t] : std::array<std::array<int, 3>, 2>{{
             {151, 211, 31584}, {197, 209, 41148}}}) {
        Int m = r * s, n = 1 + m * t, product = m * n;
        Int hi = floor_sqrt(product), limit = s * n / (s + t);
        int count = 0;
        for (Int u = (product + limit - 1) / limit; u <= hi; ++u) {
            if (product % u) continue;
            Int v = product / u;
            if (std::gcd(u, v) != 1 || (u + v - 1 - m) % t) continue;
            Row first, second;
            if (!make_row(n, u, r, s, t, first) ||
                !make_row(n, v, r, s, t, second)) continue;
            Int x = valuation(n, first), y = valuation(n, second);
            if (std::min(x, y) <= 0) throw std::runtime_error("odd example");
            ++count;
            std::cout << n << "  cuts " << u << " " << v
                      << "  v2 " << x << " " << y << '\n';
        }
        if (count != 1) throw std::runtime_error("example count");
    }
}

int main() {
    check_examples();
    std::vector<std::vector<int>> roots(bound * bound + 1);
    for (int t = 1; t <= bound * bound; ++t) roots[t] = idempotents(t);
    std::set<std::array<Int, 6>> found;
    for (int b = 2; b <= bound; ++b) {
        for (int a = 1; a < b; ++a) {
            if (std::gcd(a, b) > 1) continue;
            Int m = a * b;
            for (int t = static_cast<int>((cutoff + m - 1) / m); t <= b * b; ++t) {
                if (std::gcd(Int(t), m) > 1) continue;
                Int n = m * t + 1, product = m * n;
                Int limit = b * n / (b + t);
                Int lo = (product + limit - 1) / limit, hi = floor_sqrt(product);
                if (lo > hi) continue;
                for (int e : roots[t]) {
                    // The smaller cut is 1 or m at each full power of t.
                    Int residue = (m + (1 - m) * e - lo) % t;
                    if (residue < 0) residue += t;
                    for (Int u = lo + residue; u <= hi; u += t) {
                        if (product % u) continue;
                        Int v = product / u;
                        if (std::gcd(u, v) != 1) continue;
                        for (auto [r, s] : std::array<std::array<int, 2>, 3>{{
                                 {a, b}, {-a, b}, {a, -b}}}) {
                            Row first, second;
                            if (!make_row(n, u, r, s, t, first) ||
                                !make_row(n, v, r, s, t, second)) continue;
                            // This fixed prime suffices on every survivor of this enumeration.
                            Int x = valuation(n, first), y = valuation(n, second);
                            if (std::min(x, y) <= 0)
                                throw std::runtime_error("inspect surviving odd row");
                            if (found.insert({n, r, s, t, u, v}).second)
                                std::cout << n << "  " << r << " " << s << " " << -t
                                          << "  " << u << " " << v
                                          << "  v2 " << x << " " << y << '\n';
                        }
                    }
                }
            }
        }
    }
    std::cout << "distinct positive pairs " << found.size() << '\n';
}
