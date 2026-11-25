// Enumerate the row assigned S >= V in the endpoint factorization.
// A surviving assignment is printed for exact complementary-row verification.
#include "range_factors.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

using Wide = __int128_t;

struct Split {
    int product = 1, maximum = 0, mask = 0, count = 0;
};

std::vector<Split> splits(const Factors& f) {
    std::vector<Split> list(1 << f.size);
    for (int mask = 1; mask < int(list.size()); ++mask) {
        int j = __builtin_ctz(unsigned(mask));
        Split prev = list[mask & (mask-1)];
        int w = f.atoms[j].power;
        list[mask] = {prev.product*w, std::max(w, prev.maximum), mask, prev.count+1};
    }
    return list;
}

bool size_possible(const Split& s, int sum) {
    return sum >= 3 && 2LL*s.maximum+1 <= sum &&
           Wide(27)*s.product*s.product <= Wide(sum)*sum*sum;
}

int inverse(int a, int b) {
    int64_t r = a, q = b, u = 1, v = 0;
    while (q) {
        int64_t k = r/q, nr = r-k*q, nu = u-k*v;
        r = q; q = nr; u = v; v = nu;
    }
    if (r != 1) std::abort();
    return int((u % b + b) % b);
}

bool avoids(std::array<int,3> row, int prime) {
    // Sum of the low residues at each place. The row need not be sorted.
    while (row[0] || row[1] || row[2]) {
        int64_t total = 0;
        for (int& x : row) {
            total += x % prime;
            x /= prime;
        }
        if (total >= prime) return false;
    }
    return true;
}

struct Assignment {
    int n, u, v, s;
    std::vector<int> units, bases;
    int64_t survivors = 0;
    void visit(int j, int bins, std::array<int64_t,3> sums) {
        if (int(units.size())-j < 3-bins) return;
        if (j < int(units.size())) {
            for (int i = 0; i < std::min(3,bins+1); ++i) {
                auto next = sums;
                next[i] += units[j];
                visit(j+1, std::max(bins,i+1), next);
            }
            return;
        }
        std::array<int,3> row{};
        int64_t total = 0;
        for (int i = 0; i < 3; ++i) {
            int x = int(sums[i] % s);
            if (!x || x >= v) return;
            total += x;
            row[i] = x*u;
        }
        if (total != v) return;
        for (int p : bases) if (!avoids(row,p)) return;
        std::sort(row.begin(),row.end());
        std::cout << n << ' ' << u << ' ' << s;
        for (int x : row) std::cout << ' ' << x;
        std::cout << '\n';
        ++survivors;
    }
};

int64_t enumerate_rows(int n, const Factors& current, const Factors& pred,
                       const Split& u, const Split& s) {
    int v = n/u.product;
    Assignment a{n,u.product,v,s.product,{},{}};
    for (int j = 0; j < pred.size; ++j) {
        if (!(s.mask & (1<<j))) continue;
        auto [p,w] = pred.atoms[j];
        int m = s.product/w;
        int multiplier = int(int64_t(v%w)*inverse(m,w) % w);
        a.units.push_back(m*multiplier);
        a.bases.push_back(p);
    }
    for (int j = 0; j < current.size; ++j)
        if (u.mask & (1<<j)) a.bases.push_back(current.atoms[j].prime);
    a.visit(0,0,{});
    return a.survivors;
}

int main(int argc, char** argv) {
    // Inspect a single heavy row without filtering on its complementary row.
    if (argc == 5 && std::string(argv[1]) == "--assignment") {
        int n = std::stoi(argv[2]), scale = std::stoi(argv[3]), shift = std::stoi(argv[4]);
        if (n < 3 || n > 1000000000 || scale < 1 || shift < n/scale) return 2;
        FactorSieve sieve(n);
        auto current = sieve.factor(n), pred = sieve.factor(n-1);
        for (const auto& u : splits(current))
            for (const auto& s : splits(pred))
                if (u.product == scale && s.product == shift) {
                    enumerate_rows(n,current,pred,u,s);
                    return 0;
                }
        return 2;
    }
    if (argc != 3) return 2;
    int first = std::stoi(argv[1]), last = std::stoi(argv[2]);
    if (first < 3 || first > last || last > 1000000000) return 2;
    FactorSieve sieve(last);
    Factors pred = sieve.factor(first-1);
    int64_t survivors = 0, assignments = 0;
    for (int n = first; n <= last; ++n) {
        Factors current = sieve.factor(n);
        if (pred.size >= 3 && current.size >= 2) {
            auto scales = splits(current), shifts = splits(pred);
            for (const Split& s : shifts) {
                if (s.count < 3) continue;
                const Split& t = shifts[(shifts.size()-1)^s.mask];
                for (const Split& u : scales) {
                    int v = n/u.product;
                    if (s.product < v || !size_possible(t,u.product) || !size_possible(s,v))
                        continue;
                    ++assignments;
                    survivors += enumerate_rows(n,current,pred,u,s);
                }
            }
        }
        pred = current;
    }
    std::cerr << first << ".." << last << " size-feasible assignments " << assignments
              << " surviving rows " << survivors << '\n';
}
