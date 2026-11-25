// Exact factorizations on 2..10^9. A composite has a least divisor below 32768.
#pragma once
#include <array>
#include <cstdint>
#include <vector>

struct Atom {
    int prime, power;
};

struct Factors {
    std::array<Atom, 10> atoms{};
    int size = 0;
    const Atom* begin() const { return atoms.data(); }
    const Atom* end() const { return begin() + size; }
};

class FactorSieve {
    std::vector<uint16_t> odd_divisor;
public:
    explicit FactorSieve(int last) : odd_divisor(last / 2 + 1) {
        for (int p = 3; int64_t(p) * p <= last; p += 2) {
            if (odd_divisor[p / 2]) continue;
            for (int j = p*p; j <= last; j += 2*p)
                if (!odd_divisor[j / 2]) odd_divisor[j / 2] = p;
        }
    }
    bool prime(int n) const {
        return n == 2 || (n > 2 && n % 2 && !odd_divisor[n / 2]);
    }
    Factors factor(int n) const {
        Factors result;
        while (n > 1) {
            int p = n % 2 ? odd_divisor[n / 2] : 2;
            if (!p) {
                result.atoms[result.size++] = {n, n};
                break;
            }
            int w = 1;
            do {
                n /= p;
                w *= p;
            } while (n % p == 0);
            result.atoms[result.size++] = {p, w};
        }
        return result;
    }
};
