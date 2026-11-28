// First weights at which truncated endpoint-prime carries do not settle a pair.
#include "endpoint_splits.hpp"
#include <array>
#include <iostream>
#include <stdexcept>
#include <string>

using Row = std::array<int, 3>;
using Integer = std::int64_t;

struct Side {
    int n, scale, shift, depth, sum;
    std::vector<int> units;
    std::vector<Atom> tests;
    Row witness{};

    Side(int weight, const Factors& nf, const Factors& mf,
         const Split& u, const Split& s, int k)
        : n(weight), scale(u.product), shift(s.product), depth(k), sum(n / scale) {
        for (int j = 0; j < nf.size; ++j)
            if (u.mask & (1 << j)) tests.push_back(nf.atoms[j]);
        for (int j = 0; j < mf.size; ++j) {
            if (!(s.mask & (1 << j))) continue;
            auto atom = mf.atoms[j];
            int other = shift / atom.power;
            int a = Integer(sum % atom.power) * inverse(other, atom.power) % atom.power;
            units.push_back(a * other);
            tests.push_back(atom);
        }
    }

    bool passes(const Row& row) const {
        for (auto [p, q] : tests) {
            Integer mod = q;
            for (int j = 0; j < depth && mod <= n / p; ++j) {
                mod *= p;
                if (row[0] % mod + row[1] % mod + row[2] % mod != n % mod)
                    return false;
            }
        }
        return true;
    }

    bool lifts(const std::array<Integer, 3>& bins) {
        Row low;
        for (int j = 0; j < 3; ++j)
            low[j] = (bins[j] + shift - 1) % shift + 1;
        Integer spare = Integer(sum) - low[0] - low[1] - low[2];
        if (spare < 0) return false;
        if (spare % shift) throw std::logic_error("incompatible row sum");
        for (Integer i = 0; i <= spare; i += shift)
            for (Integer j = 0; j <= spare - i; j += shift) {
                Row row{int(scale * (low[0] + i)),
                        int(scale * (low[1] + j)),
                        int(scale * (low[2] + spare - i - j))};
                if (passes(row)) {
                    std::sort(row.begin(), row.end());
                    witness = row;
                    return true;
                }
            }
        return false;
    }

    bool visit(int next, int used, std::array<Integer, 3> bins) {
        if (shift >= sum && int(units.size()) - next < 3 - used) return false;
        if (next == int(units.size())) return lifts(bins);
        for (int i = 0; i < std::min(used + 1, 3); ++i) {
            auto updated = bins;
            updated[i] += units[next];
            if (visit(next + 1, std::max(used, i + 1), updated)) return true;
        }
        return false;
    }
};

bool find_pair(int n, const Factors& nf, const Factors& mf, int k) {
    if (nf.size < 2 || mf.size < 3) return false;
    auto weights = splits(nf), shifts = splits(mf);
    for (const auto& s : shifts) {
        if (s.count < 3) continue;
        const Split& t = shifts[(shifts.size() - 1) ^ s.mask];
        for (const auto& u : weights) {
            const Split& v = weights[(weights.size() - 1) ^ u.mask];
            if (s.product < v.product || !size_possible(s, v.product) ||
                !size_possible(t, u.product)) continue;
            Side a(n, nf, mf, u, s, k);
            if (!a.visit(0, 0, {})) continue;
            Side b(n, nf, mf, v, t, k);
            if (!b.visit(0, 0, {})) continue;
            std::cout << k << " " << n;
            for (int x : a.witness) std::cout << ' ' << x;
            for (int x : b.witness) std::cout << ' ' << x;
            std::cout << std::endl;
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc != 4) throw std::invalid_argument("expected first last initial_depth");
    int first = std::stoi(argv[1]), last = std::stoi(argv[2]), k = std::stoi(argv[3]);
    if (first < 3 || last < first || last > 1000000000 || k < 0 || k > 30)
        throw std::invalid_argument("range or depth");
    FactorSieve sieve(last);
    auto pred = sieve.factor(first - 1);
    for (int n = first; n <= last; ++n) {
        auto current = sieve.factor(n);
        while (find_pair(n, current, pred, k)) {
            if (++k > 30) throw std::runtime_error("inspect full counterexample");
        }
        pred = current;
    }
    std::cout << "range " << first << ' ' << last << " sufficient depth " << k << '\n';
}
