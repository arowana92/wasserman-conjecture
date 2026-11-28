// Independent calculation in base 5. No arbitrary-precision library is needed.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using U = std::uint64_t;
constexpr int width = 24;
constexpr U block = 59604644775390625ULL; // 5^24
constexpr U ternary = 1853020188851841ULL; // 3^32

U multiply(U a, U b, U modulus) {
    return static_cast<__uint128_t>(a) * b % modulus;
}

// Work modulo block^(words.size()). Each block is 1 modulo 2 and modulo 3.
template<int divisor>
void divide_modulo(std::vector<U>& words) {
    U sum = 0;
    for (U value : words) sum += value;
    U carry = (divisor - sum % divisor) % divisor;
    for (auto i = words.rbegin(); i != words.rend(); ++i) {
        U value = carry * block + *i;
        *i = value / divisor;
        carry = value % divisor;
    }
    assert(!carry);
}

bool fits(U part, U total, U base) {
    while (part) {
        if (part % base > total % base) return false;
        part /= base;
        total /= base;
    }
    return true;
}

U digit_slice(const std::vector<unsigned char>& digits, int low, int high, U modulus) {
    U value = 0;
    for (int i = high-1; i >= low; --i)
        value = (5*value + digits[i]) % modulus;
    return value;
}

int main(int argc, char** argv) {
    int limit = argc == 2 ? std::atoi(argv[1]) : 1024;
    if (limit < 1 || limit > 1024) return 2;
    int blocks = (limit + width + width-1)/width;
    std::vector<U> inverse_r(blocks);
    inverse_r[0] = 1;
    std::vector<unsigned char> digits(width*blocks);
    std::vector<U> r5(limit+1), d3(limit+1, 1);
    for (int a = 1; a <= limit; ++a) d3[a] = d3[a-1]*5 % ternary;
    U r2 = 1, r_mod_block = 1;
    std::array<U, 3> counts{};
    int highest = 0;
    for (int c = 1; c <= limit; ++c) {
        divide_modulo<3>(inverse_r);
        r2 *= 3; // unsigned wrap computes residues modulo 2^64
        r_mod_block = r_mod_block * 3 % block;
        r5[0] = r2;
        for (int a = 1; a <= limit; ++a) r5[a] = 5*r5[a-1];
        auto inverse = inverse_r;
        U q3 = 1, w5 = r_mod_block;
        for (int b = 1; b <= limit; ++b) {
            divide_modulo<2>(inverse);
            q3 = 2*q3 % ternary;
            w5 = 2*w5 % block;
            int pos = 0;
            for (U chunk : inverse) {
                for (int i = 0; i < width; ++i) {
                    digits[pos++] = static_cast<unsigned char>(chunk%5);
                    chunk /= 5;
                }
            }
            U u2 = 0;
            for (int a = 1; a <= limit; ++a) {
                U digit = digits[a-1];
                u2 += digit*r5[a-1];
                if (u2 & ~r5[a]) continue;
                ++counts[0];
                U x3 = digit_slice(digits, 0, a, ternary);
                if (!fits(multiply(x3, q3, ternary), multiply(d3[a], q3, ternary), 3))
                    continue;
                ++counts[1];
                U next = digit_slice(digits, a, a+width, block);
                U quotient = (block - multiply(next, w5, block)) % block;
                if (fits(quotient, w5, 5)) {
                    ++counts[2];
                    std::cout << "passing digits " << a << ' ' << b << ' ' << c << '\n';
                } else {
                    U modulus = 5;
                    int depth = 1;
                    while (quotient % modulus <= w5 % modulus) {
                        modulus *= 5;
                        ++depth;
                    }
                    highest = std::max(highest, a+depth);
                    std::cout << "base5 " << a << ' ' << b << ' ' << c << ' '
                              << depth << '\n';
                }
            }
        }
    }
    std::cout << "maximum base-5 position " << highest << '\n';
    std::cout << U(limit)*limit*limit;
    for (auto n : counts) std::cout << ' ' << n;
    std::cout << std::endl;
    if (limit == 1024 && (counts[0] != 3182 || counts[1] != 7 || highest != 9)) return 1;
    return counts[2] != 0;
}
