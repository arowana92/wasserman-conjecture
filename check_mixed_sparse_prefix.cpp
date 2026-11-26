// Compile with a C++17 compiler supporting unsigned 128-bit integers.
// Certify the high binary digits of possible sums of six powers of p.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <utility>

using Word = __uint128_t;

int high_bit(Word x) {
    const auto upper = static_cast<std::uint64_t>(x >> 64);
    return upper ? 127 - __builtin_clzll(upper)
                 : 63 - __builtin_clzll(static_cast<std::uint64_t>(x));
}

int popcount(Word x) {
    return __builtin_popcountll(static_cast<std::uint64_t>(x >> 64))
         + __builtin_popcountll(static_cast<std::uint64_t>(x));
}

// Least integer >= x with at most six ones.
Word next_sparse(Word x) {
    if (popcount(x) <= 6) return x;
    Word tail = x;
    Word bit = 0;
    for (int i = 0; i < 6; ++i) {
        bit = Word(1) << high_bit(tail);
        tail -= bit;
    }
    return x - tail + bit;
}

struct Interval {
    Word lo;
    Word hi;
};

Interval divide_power_two(Interval v, int shift) {
    const Word mask = (Word(1) << shift) - 1;
    return {v.lo >> shift, (v.hi + mask) >> shift};
}

class PrefixSearch {
    std::array<Interval, 105> powers{};

    bool possible(Interval prefix, int last, unsigned left) const {
        const Word target = next_sparse(prefix.lo);
        if (target <= prefix.hi + left) return true;  // Remaining powers in [0,1].
        if (!left || target > prefix.hi + left * powers[last].hi) return false;
        const Word needed = (target-prefix.hi+left-1)/left;
        // Search in descending order so a whole interval of small summands
        // is discarded as soon as even their common upper bound is too low.
        for (int j = last; j > 0 && powers[j].hi >= needed; --j) {
            Interval next{prefix.lo+powers[j].lo, prefix.hi+powers[j].hi};
            if (possible(next, j, left-1)) return true;
        }
        return false;
    }

public:
    bool test(int p, Interval leading) {
        // Keep the leading power and at most 100 of its highest bits.
        int shift = std::max(0, high_bit(leading.hi)-99);
        Interval term = divide_power_two(leading, shift);
        std::array<Interval, 105> descending{};
        int count = 0;
        while (term.hi > 1) {
            descending[count++] = term;
            term = {term.lo / p, (term.hi+p-1) / p};
        }
        powers[0] = {0, 1};
        const int top = count;
        for (int i = 0; i < count; ++i) powers[top-i] = descending[i];
        return possible(powers[top], top, 5);
    }
};

bool prime(int p) {
    for (int q = 2; q*q <= p; ++q) if (p%q == 0) return false;
    return true;
}

void positive_controls() {
    for (Word n = 0; n <= 32768; ++n) {
        Word expected = n;
        while (popcount(expected) > 6) ++expected;
        assert(next_sparse(n) == expected);
    }
    // Make the leading power of each actual positive example a 116-bit
    // mantissa by shifting all six summands by the same binary exponent.
    const std::array<std::pair<int, std::uint64_t>, 4> examples{{
        {3, 10737942553ULL}, {5, 8594141251ULL}, {7, 571539529ULL}, {251, 1507ULL}
    }};
    PrefixSearch sieve;
    for (const auto& example : examples) {
        int p = example.first;
        auto n = example.second-1;
        int digit_sum = 0;
        Word leading = 1;
        for (auto value = n; value; value /= p) {
            digit_sum += value % p;
            if (value >= static_cast<unsigned>(p)) leading *= p;
        }
        assert(digit_sum == 6 && popcount(n) <= 6);
        leading <<= 115-high_bit(leading);
        assert(sieve.test(p, {leading, leading}));
    }
}

int main(int argc, char** argv) {
    const long long bits = argc == 2 ? std::atoll(argv[1]) : 1000000;
    if (argc > 2 || bits < 128 || bits > 100000000) return 2;
    positive_controls();
    PrefixSearch sieve;
    long long checked = 0, survivors = 0;
    for (int p = 3; p < 256; p += 2) {
        if (!prime(p)) continue;
        Interval mantissa{1, 1};
        long long scale = 0;
        for (long long e = 1; ; ++e) {
            mantissa.lo *= p;
            mantissa.hi *= p;
            const int cut = std::max(0, high_bit(mantissa.hi)-115);
            mantissa = divide_power_two(mantissa, cut);
            scale += cut;
            if (scale+high_bit(mantissa.lo) >= bits) break;
            if (scale+high_bit(mantissa.hi) < 120) continue;
            ++checked;
            if (sieve.test(p, mantissa)) {
                ++survivors;
                std::printf("Unexcluded leading power %d^%lld\n", p, e);
            }
        }
        std::fprintf(stderr, "Through prime %d: %lld prefixes, %lld survivors\n",
                     p, checked, survivors);
    }
    std::printf("%lld prefixes, %lld survivors, ceiling 2^%lld.\n", checked, survivors, bits);
    return survivors != 0;
}
