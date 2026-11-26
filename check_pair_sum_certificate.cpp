// Exhaust every positive binary submask in the two four-prime certificates.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

int main() {
    // N = 2^80 + offset. The offset fits in 64 bits.
    constexpr std::uint64_t offset = 47620826697754176ULL;
    constexpr std::array<std::uint64_t, 2> primes = {134217689, 268947953};
    std::array<std::vector<std::uint64_t>, 2> zero, one;
    std::uint64_t visited = 0;
    for (auto part = offset; part > 0; part = (part - 1) & offset) {
        ++visited;
        for (unsigned i = 0; i < primes.size(); ++i) {
            const auto residue = part % primes[i];
            if (residue == 0) zero[i].push_back(part);
            if (residue == 1) one[i].push_back(part);
        }
    }
    if (visited != (1ULL << 29) - 1) return 1;
    for (unsigned i = 0; i < primes.size(); ++i) {
        if (zero[i].size() != 4 || one[i].size() != i) return 1;
        for (auto a : zero[i]) {
            for (auto b : zero[i]) if ((a & b) == 0) return 1;
            for (auto b : one[i]) if ((a & b) == 0) return 1;
        }
    }
    std::cout << visited << " positive submasks checked against both primes.\n";
}
