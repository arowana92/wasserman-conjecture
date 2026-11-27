// Enumerate the finite range for two single prime-power summands.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using I = std::int64_t;
using Wide = __int128_t;
constexpr I checked_interval = 1000000000;
constexpr I height_limit = 2000000000000; // Above 131072^(5/2)/sqrt(27) + 1.

Wide cube(I n) { return Wide(n) * n * n; }
Wide fifth(I n) { return cube(n) * n * n; }

bool eligible(I n, I m, I p, I high, I q, I low) {
    I pred = n - 1;
    Wide target = Wide(27) * pred * pred;
    bool split_bound = q > low && (high > low || p != 2);
    if (q <= low && fifth(m * q) <= target) return false;
    if (q > 2 * low && high == low && p == 2 &&
        cube(low) * m * m * m * m < 2 * target) return false;

    I mask = 0;
    // For an odd base, store the nonzero digit steps and their maximum digits.
    std::array<I, 32> steps{}, digits{}, limits{};
    int used = 0;
    if (p == 2) {
        mask = m & (n / q);
    } else {
        I other = n / q;
        for (I val = m, step = 1; val; val /= p, other /= p, step *= p) {
            I digit = std::min(val % p, other % p);
            if (digit) {
                limits[used] = digits[used] = digit;
                steps[used++] = step;
                mask += digit * step;
            }
        }
    }
    for (I k = mask; k > 0;) {
        I edge = m + (low - 1) * k;
        if (split_bound && cube(low * k) * edge * edge < target) break;
        if (q <= low && fifth(k * q) <= target) break;
        I z = k * q, u = std::gcd(n, z);
        if (z >= 3 && n - z >= 6 && u >= 3) {
            I odd = std::gcd(n / low, k);
            I bound = k * edge;
            if (!split_bound || cube(low) * odd * bound * bound >= target) {
                I s = std::gcd(pred, z) * std::gcd(pred, z - 1);
                I t = pred / s;
                if (cube(u) >= Wide(27) * t * t) return true;
            }
        }
        if (p == 2) {
            k = (k - 1) & mask;
        } else {
            int i = 0;
            while (!digits[i]) {
                digits[i] = limits[i];
                k += digits[i] * steps[i];
                ++i;
            }
            --digits[i];
            k -= steps[i];
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    I cap = std::atoll(argv[1]);
    std::string mode = argv[2];
    if (cap < 2 || cap > 131072 ||
        (mode != "binary-ordered" && mode != "binary-above" && mode != "odd")) return 2;
    std::vector<bool> composite(cap + 1);
    I families = 0, weights = 0, feasible = 0, maximum = 0;
    for (I p = 2; p <= cap; ++p) {
        if (composite[p]) continue;
        for (I mult = p*p; mult <= cap; mult += p) composite[mult] = true;
        if ((p == 2) == (mode == "odd")) continue;
        std::vector<I> powers{1};
        while (powers.back() <= height_limit / p) powers.push_back(powers.back() * p);
        for (I low : powers) {
            if (low == 1 || low > cap) continue;
            for (I m = 1; m <= cap / low; ++m) {
                ++families;
                Wide bound = fifth(low * m);
                for (I high : powers) {
                    if (high < low) continue;
                    if (Wide(27) * high * high > bound) break;
                    for (I q : powers) {
                        if (mode == "binary-ordered" && q > high) break;
                        if (mode == "binary-above" && q <= high) continue;
                        if (q == 1 && m % p > 1) continue;
                        I n = high + m * q + low;
                        Wide target = Wide(27) * (n - 1) * (n - 1);
                        if (target > bound) break;
                        if (n <= checked_interval) continue;
                        ++weights;
                        maximum = std::max(maximum, n);
                        if (eligible(n, m, p, high, q, low)) {
                            ++feasible;
                            std::cout << p << ' ' << m << ' ' << high << ' ' << q
                                      << ' ' << low << ' ' << n << '\n';
                        }
                    }
                }
            }
        }
    }
    std::cerr << families << " families, " << weights << " parameter choices, "
              << feasible << " with an eligible part, maximum weight " << maximum << '\n';
}
