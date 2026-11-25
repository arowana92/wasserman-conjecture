// Search CRT rows with two specified prime factors and a prime scale.
// The scale and the two factors are prime. The last factor may be any integer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Integer = std::int64_t;
using Counts = std::array<Integer, 6>;

Integer inverse(Integer value, Integer modulus) {
    Integer b = modulus, x = 1, y = 0;
    while (b) {
        Integer q = value / b;
        value -= q*b;
        std::swap(value,b);
        x -= q*y;
        std::swap(x,y);
    }
    if (value != 1) throw std::logic_error("nonunit in CRT progression");
    return (x%modulus+modulus)%modulus;
}

bool carry_free(Integer a, Integer b, Integer c, int prime) {
    while (a || b || c) {
        if (a%prime+b%prime+c%prime >= prime) return false;
        a /= prime;
        b /= prime;
        c /= prime;
    }
    return true;
}

std::vector<int> odd_primes(int limit) {
    std::vector<int> result;
    for (int n=3; n<=limit; n+=2) {
        bool prime = true;
        for (int p : result) {
            if (p*p>n) break;
            if (n%p==0) { prime=false; break; }
        }
        if (prime) result.push_back(n);
    }
    return result;
}

Counts search(int limit, int first_scale, int last_scale, int periods) {
    // All products, including each weight, fit a signed 64-bit integer.
    if (limit<3 || limit>3659 || first_scale<3 || last_scale<first_scale ||
        last_scale>800 || periods<1 || periods>12)
        throw std::invalid_argument("bounds exceed the certified arithmetic range");
    Counts counts{};
    auto primes = odd_primes(limit);
    for (int u : primes) {
        if (u<first_scale || u>last_scale) continue;
        for (int p : primes) {
            if (p<=u) continue;
            for (int q : primes) {
                if (q<=p) continue;
                Integer pq = Integer(p)*q, period = u*pq;
                for (int a=1; u*a<p; ++a) {
                    Integer left = Integer(u)*a*q;
                    for (int b=1; left+Integer(u)*b*p<pq; ++b) {
                        Integer middle = Integer(u)*b*p, right = pq-left-middle;
                        Integer start = period-inverse(right,period);
                        ++counts[0];
                        int skip = start<=q ? 1 : 0;
                        counts[2] += periods-skip;
                        // These three inequalities depend only on the progression.
                        if ((a*q*(start%u))%u+(b*p*(start%u))%u>=u ||
                            (u*b*(start%p))%p > (q*(start%p))%p ||
                            (u*a*(start%q))%q > (p*(start%q))%q) continue;
                        ++counts[1];
                        for (int j=skip; j<periods; ++j) {
                            Integer r = start+period*j;
                            Integer n = pq*r+1;
                            Integer A = left*r, B = middle*r, C = right*r+1;
                            if ((n/u)%u==0 || !carry_free(A/u,B/u,C/u,u)) continue;
                            ++counts[3];
                            if (!carry_free(A,B,C,p)) continue;
                            ++counts[4];
                            if (!carry_free(A,B,C,q)) continue;
                            ++counts[5];
                            std::cout << "candidate " << u << ' ' << p << ' ' << q << ' '
                                      << r << ' ' << A << ' ' << B << ' ' << C << '\n';
                        }
                    }
                }
            }
        }
    }
    return counts;
}

int main(int argc, char** argv) {
    if (argc!=1 && argc!=5) {
        std::cerr << "usage: check_three_factor_progressions [limit first_u last_u periods]\n";
        return 2;
    }
    Counts counts = argc==1 ? search(3659,103,800,12)
        : search(std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]));
    std::cout << "progressions, local survivors, parameters, scale survivors, "
                 "first-base survivors, final survivors\n";
    for (auto n : counts) std::cout << n << ' ';
    std::cout << '\n';
}
