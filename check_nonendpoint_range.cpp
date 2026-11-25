// Exclude selections without (1,1,N-2) using two-row weighted products.
#include "range_factors.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <string>

constexpr int max_height = 30;
constexpr int ring_size = max_height + 1;
int triangle(int t) { return t*(t+1)/2; }

// Lower and upper integer bounds on 8 log_2(n). Entries fit in 64 bits
// since the table mantissas are at most 255.
struct LogBounds {
    std::array<int,256> mantissa{};
    LogBounds() {
        for (uint64_t m = 128; m < 256; ++m) {
            uint64_t sq = m*m, fourth = sq*sq, eighth = fourth*fourth;
            mantissa[m] = 63-__builtin_clzll(eighth)-56;
        }
    }
    int low(int n) const {
        int e = 31-__builtin_clz(unsigned(n));
        int m = e >= 7 ? n >> (e-7) : n << (7-e);
        return 8*e + mantissa[m];
    }
    int high(int n) const { return low(n) + ((n & (n-1)) ? 2 : 0); }
} logs8;

int bounded_part(int m) {
    int cost = 0;
    for (int i = 1; i <= m; ++i) cost += i*logs8.high(i);
    return cost;
}

int row_bound(int log_n, int h, int min_part) {
    int t = (h+2)/3;
    int degree = 0;
    for (int i = 0; i < 3; ++i) degree += triangle((h+i)/3);
    int bound = degree*(log_n-logs8.low(3));
    for (int m = min_part; m < t; ++m) {
        int u = (h-m+1)/2, v = h-m-u;
        int cost = bounded_part(m);
        bound = std::max(bound,cost+(triangle(u)+triangle(v))*(log_n-8));
        for (int y = std::max(m,11-m); y < u; ++y) {
            int rest = h-m-y;
            bound = std::max(bound,cost+bounded_part(y)+triangle(rest)*log_n);
        }
    }
    return bound;
}

bool carries(int n, int a, int b, int p) {
    for (int64_t w = p; w <= n; w *= p)
        if (a % w + b % w > n % w) return true;
    return false;
}

// Lower bound for the contribution at one p-power layer from three cutoffs whose
// differences from the respective residues sum to at least d.
int layer_cost(int d, int w) {
    int result = 0;
    for (; d > 0; d -= 3*w) result += d;
    return result;
}

// With exact_sum false, small_sum is a lower bound for a+b. With a zero
// first part and an exact sum, only binomial(N, small_sum) is used.
bool excludes(int n, int gap, int small_sum, int a, int height,
              const std::array<Factors,ring_size>& factors, bool exact_sum = true) {
    static std::map<std::array<int,3>,int> cache;
    int log_n = logs8.high(n);
    int common = exact_sum ? std::gcd(n,small_sum) : n;
    if (a) common = std::gcd(common,a);
    int minimum = common <= 2 ? 3 : 1;
    std::array<int,3> key{log_n,height,minimum};
    auto found = cache.find(key);
    if (found == cache.end()) {
        int upper = row_bound(log_n,height,3) + row_bound(log_n,height,minimum);
        found = cache.emplace(key,upper).first;
    }
    int bound = found->second, lower = 0;
    for (int p : {2,3,5,7,11,13,17,19,23,29}) {
        if (p > height) break;
        int64_t first_above_gap = p;
        while (first_above_gap <= gap) first_above_gap *= p;
        bool forced = n % first_above_gap < small_sum ||
                      (exact_sum && carries(n,a,small_sum-a,p));
        for (int w = p; w <= height; w *= p) {
            int h = n%w;
            int largest_sum = h + w*((3*w-3-h)/w);
            int baseline = layer_cost(height-largest_sum,w);
            int other = forced ? layer_cost(height-h,w) : baseline;
            lower += (baseline+other)*logs8.low(p);
        }
    }
    for (int h = 0; h < height; ++h) {
        int part = 1;
        for (auto [p,w] : factors[(n-h)%ring_size]) {
            if (w <= height) continue;
            if (!(w > gap && h < w && h < small_sum) &&
                !(exact_sum && carries(n,a,small_sum-a,p))) continue;
            int allowed = w;
            for (int k = height; k >= p; k /= p) allowed /= p;
            part *= allowed;
        }
        lower += (height-h)*logs8.low(part);
        if (lower > bound) return true;
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int first = std::stoi(argv[1]), last = std::stoi(argv[2]);
    if (first < 3*max_height || first > last || last > 1000000000) return 2;
    FactorSieve sieve(last);
    std::array<Factors,ring_size> factors{};
    for (int n = first-max_height; n < first; ++n)
        factors[n%ring_size] = sieve.factor(n);
    int prime = first;
    while (!sieve.prime(prime)) --prime;
    int64_t short_gap = 0, whole = 0, by_sum = 0, by_parts = 0, left = 0;
    for (int n = first; n <= last; ++n) {
        factors[n%ring_size] = sieve.factor(n);
        if (sieve.prime(n)) prime = n;
        int gap = n-prime;
        if (gap < 11) { ++short_gap; continue; }
        if (excludes(n,gap,11,0,9,factors,false) || excludes(n,gap,11,0,11,factors,false)) {
            ++whole;
            continue;
        }
        for (int sum = 11; sum <= gap; ++sum) {
            int height = std::min(sum,max_height);
            if (excludes(n,gap,sum,0,height,factors)) {
                by_sum += sum/2;
                continue;
            }
            for (int a = 1; a <= sum/2; ++a) {
                bool ok = false;
                for (int h : {6,9,11,height,18,30})
                    if (excludes(n,gap,sum,a,h,factors)) { ok = true; break; }
                if (ok) ++by_parts;
                else { ++left; std::cout << n << ' ' << a << ' ' << sum-a << '\n'; }
            }
        }
    }
    std::cerr << first << ".." << last << " short-gap weights " << short_gap
              << " whole-weight product exclusions " << whole << " rows by sum " << by_sum
              << " rows by parts " << by_parts << " remaining rows " << left << '\n';
}
