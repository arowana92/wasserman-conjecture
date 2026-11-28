// Enumerate cut pairs by the smaller term in their reduced ratio.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

using Int = std::int64_t;

bool avoids(Int n, Int a, Int p) {
    while (a) {
        if (a%p > n%p) return false;
        a /= p;
        n /= p;
    }
    return true;
}

bool shares_prime(Int n, Int a, Int b, Int p) {
    return n%p <= 1 && !avoids(n,a,p) && !avoids(n,b,p);
}

bool passes_factors(Int n, Int a, Int b, Int factor) {
    for (Int p = 2; p*p <= factor; p = p == 2 ? 3 : p+2) {
        if (factor%p) continue;
        if (shares_prime(n,a,b,p)) return false;
        do { factor /= p; } while (factor%p == 0);
    }
    return factor == 1 || !shares_prime(n,a,b,factor);
}

std::vector<int> least_factors(int limit) {
    std::vector<int> least(limit+1);
    for (int p = 2; p <= limit; ++p) {
        if (least[p]) continue;
        for (int j = p; j <= limit; j += p) {
            if (!least[j]) least[j] = p;
        }
    }
    return least;
}

void divisors(int value, const std::vector<int>& least, std::vector<int>& result) {
    result.assign(1,1);
    while (value > 1) {
        int p = least[value], power = 1;
        auto initial = result.size();
        do {
            value /= p;
            power *= p;
            for (std::size_t i = 0; i < initial; ++i) result.push_back(result[i]*power);
        } while (value%p == 0);
    }
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) return 2;
    int limit = std::atoi(argv[1]);
    Int cutoff = std::atoll(argv[2]);
    bool audit = argc == 4 && std::string_view(argv[3]) == "--audit";
    if (limit < 2 || limit > 4096 || cutoff < 0 || cutoff > 1000000000 ||
        (argc == 4 && !audit)) return 2;
    auto least = least_factors(limit*limit/4);
    std::vector<int> options;
    Int systems = 0, odd = 0, passed = 0, largest = 0;
    for (int x = 2; x <= limit; ++x) {
        for (int m = 1; m < x; ++m) {
            Int u = Int(m)*(x-m), top = m*u;
            int alpha = std::gcd(x,m);
            for (int k = (u-1)/x; k >= 1; --k) {
                Int L = u-Int(k)*x;
                // h*N = k*x*y, while L*y = top-h. Lowering k only lowers
                // this upper bound on h and increases L.
                Int hmax = std::min<Int>(k*m, (Int(k)*x*top)/(L*(cutoff+1)));
                if (hmax < alpha || top <= L*x) break;

                auto check = [&](int h) {
                    Int y = (top-h)/L;
                    if (y < x || std::gcd<Int>(x,y) != 1 || Int(k)*x*y%h) return;
                    Int g = k*m/h, n = Int(k)*x*y/h, a = g*x, b = g*y;
                    if (a < 2 || b >= n || n <= cutoff) return;
                    ++systems;
                    largest = std::max(largest,n);
                    if ((a & (n-a)) && (b & (n-b))) return;
                    ++odd;
                    if (audit) std::cout << "ODD " << x << ' ' << y << ' ' << n
                                         << ' ' << a << ' ' << b << '\n';
                    // h*N=k*x*y and h*(N-1)=m*(x-m)*(y-m) cover every prime.
                    for (Int f : {Int(k),Int(x),y,Int(m),Int(x-m),y-m}) {
                        if (!passes_factors(n,a,b,f)) return;
                    }
                    ++passed;
                    std::cout << "PASS " << x << ' ' << y << ' ' << n << ' '
                              << a << ' ' << b << '\n';
                };

                Int first = top%L;
                if (!first) first = L;
                if (first > hmax) continue;
                if ((hmax-first)/L < 32) {
                    for (Int h = first; h <= hmax; h += L) {
                        if (h%alpha == 0 && k*m%h == 0) check(int(h));
                    }
                } else {
                    // alpha divides h, and h divides k*m.
                    divisors(k*(m/alpha),least,options);
                    for (int d : options) {
                        int h = alpha*d;
                        if (h <= hmax && (top-h)%L == 0) check(h);
                    }
                }
            }
        }
    }
    std::cout << "Systems " << systems << ", passing parity " << odd
              << ", passing all bases " << passed << ", largest weight " << largest << '\n';
}
