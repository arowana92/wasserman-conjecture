// Enumerate q first, using (m-1)^2 = 1 mod q.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using Int = std::int64_t;

Int inverse(Int a, Int b) {
    Int modulus = b, r = 1, s = 0;
    while (b) {
        Int k = a/b, next = a-k*b, coefficient = r-k*s;
        a = b;
        b = next;
        r = s;
        s = coefficient;
    }
    return (r%modulus+modulus)%modulus;
}

std::vector<Int> roots_of_one(int q, const std::vector<int>& least) {
    std::vector<Int> roots{0};
    Int modulus = 1;
    while (q > 1) {
        int p = least[q], power = 1;
        do { q /= p; power *= p; } while (q%p == 0);
        std::vector<int> local{1};
        if (power > 2) local.push_back(power-1);
        if (p == 2 && power >= 8) {
            local.push_back(power/2-1);
            local.push_back(power/2+1);
        }
        std::vector<Int> next;
        Int inv = inverse(modulus,power);
        for (Int a : roots) for (Int b : local) {
            Int correction = ((b-a)*inv%power+power)%power;
            next.push_back(a+modulus*correction);
        }
        roots.swap(next);
        modulus *= power;
    }
    return roots;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int limit = std::atoi(argv[1]);
    if (limit < 3 || limit > 1000000) return 2;
    std::vector<int> least(limit+1);
    for (int p = 2; p <= limit; ++p) {
        if (least[p]) continue;
        for (int j = p; j <= limit; j += p) {
            if (!least[j]) least[j] = p;
        }
    }
    Int count = 0, odd = 0, largest = 0;
    for (int q = 1; q+3 <= limit; ++q) {
        for (Int r : roots_of_one(q,least)) {
            Int start = q+3+((r-2)%q+q)%q;
            for (Int m = start; m <= limit; m += q) {
                Int n = 1+m*(m-2)/q;
                Int disc = 1+4*m*n;
                if (disc%16 != 1 && disc%16 != 9) continue;
                Int root = std::sqrt(static_cast<long double>(disc));
                while ((root+1)*(root+1) <= disc) ++root;
                while (root*root > disc) --root;
                if (root*root != disc) continue;
                Int a = (root-1)/2, b = a+1;
                if (b >= n) continue;
                ++count;
                largest = std::max(largest,n);
                if ((a & (n-a)) && (b & (n-b))) continue;
                ++odd;
                std::cout << m << ' ' << n << ' ' << a << ' ' << b << '\n';
            }
        }
    }
    std::cout << "Systems " << count << " parity " << odd << " largest " << largest << '\n';
}
