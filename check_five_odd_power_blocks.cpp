// Verify the first 120 binary columns, pruning a prefix as soon as at most
// two exponents remain. Output any longer branches for arbitrary-precision checking.
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using Word = __uint128_t;
constexpr int cutoff = 120;

Word power(Word x, unsigned exponent) {
    Word result = 1;
    while (exponent) {
        if (exponent & 1U) result *= x;
        x *= x;
        exponent >>= 1;
    }
    return result;
}

std::string decimal(Word value) {
    if (value < 10) return std::string(1, char('0' + value));
    return decimal(value / 10) + char('0' + value % 10);
}

struct Search {
    std::array<std::array<std::vector<unsigned>, 2>, cutoff> children;
    std::array<std::uint64_t, cutoff + 1> active{};

    void visit(int bits, Word residue, const std::vector<unsigned>& exponents) {
        ++active[bits];
        if (bits == cutoff) {
            std::cout << "suffix " << decimal(residue);
            for (unsigned e : exponents) std::cout << ' ' << e;
            std::cout << '\n';
            return;
        }
        auto& low = children[bits][0];
        auto& high = children[bits][1];
        low.clear();
        high.clear();
        Word bit = Word(1) << bits;
        for (unsigned e : exponents) {
            Word q = power(residue, e);
            // Raising an odd input to an odd exponent flips the new bit
            // when the new input bit flips. The digit of 2*q stays fixed.
            if (!(2*q & bit)) {
                low.push_back(e);
                high.push_back(e);
            } else if (3*q & bit) {
                high.push_back(e);
            } else {
                low.push_back(e);
            }
        }
        if (low.size() >= 3) visit(bits + 1, residue, low);
        if (high.size() >= 3) visit(bits + 1, residue + bit, high);
    }
};

// Independent traversal. Index the possible bases by q = 3^u with u odd.
// Extend u one bit at a time and carry each entire modular power as data.
// The test compares the residues of 3*q^e and 5*q^e at a factorial layer.
struct LogSearch {
    using Entry = std::pair<unsigned, Word>;
    std::array<std::vector<Word>, cutoff> factors;
    std::array<Word, cutoff> steps{};
    std::array<std::array<std::vector<Entry>, 2>, cutoff> children;
    std::array<std::uint64_t, cutoff + 1> active{};

    explicit LogSearch(int count) {
        Word step = 9;
        for (int j = 3; j < cutoff; ++j) {
            steps[j] = step;
            Word f = step;
            for (int i = 0; i < count; ++i) {
                factors[j].push_back(f);
                f *= step*step;
            }
            step *= step;
        }
    }

    void visit(int bits, Word base, const std::vector<Entry>& powers) {
        ++active[bits];
        if (bits == cutoff) {
            std::cout << "suffix " << decimal(base % (Word(1) << cutoff));
            for (const auto& entry : powers) std::cout << ' ' << entry.first;
            std::cout << '\n';
            return;
        }
        Word mask = (Word(2) << bits)-1;
        for (int branch = 0; branch < 2; ++branch) {
            auto& next = children[bits][branch];
            next.clear();
            for (auto [e, value] : powers) {
                if (branch) value *= factors[bits][(e-1)/2];
                if ((3*value & mask) <= (5*value & mask)) next.emplace_back(e,value);
            }
            if (next.size() >= 3)
                visit(bits+1, branch ? base*steps[bits] : base, next);
        }
    }
};

int main(int argc, char** argv) {
    const int count = argc > 1 ? std::stoi(argv[1]) : 1000;
    if (count < 3 || count > 1000) return 2;
    std::array<std::uint64_t, cutoff+1> levels;
    if (argc == 3 && std::string(argv[2]) == "--logs") {
        LogSearch search(count);
        std::vector<LogSearch::Entry> powers;
        Word v = 3;
        for (int i = 0; i < count; ++i) {
            powers.emplace_back(2*i+1,v);
            v *= 9;
        }
        search.visit(3,3,powers);
        levels = search.active;
    } else if (argc <= 2) {
        Search search;
        std::vector<unsigned> exponents;
        for (int i = 0; i < count; ++i) exponents.push_back(2*i+1);
        search.visit(3,3,exponents);
        levels = search.active;
    } else {
        return 2;
    }
    for (int i = 3; i <= cutoff; ++i)
        std::cout << "level " << i << ' ' << levels[i] << '\n';
}
