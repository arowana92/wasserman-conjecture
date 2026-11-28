#pragma once
#include "range_factors.hpp"
#include <algorithm>
#include <cstdlib>

using Wide = __int128_t;

struct Split {
    int product = 1, maximum = 0, mask = 0, count = 0;
};

inline std::vector<Split> splits(const Factors& f) {
    std::vector<Split> list(1 << f.size);
    for (int mask = 1; mask < int(list.size()); ++mask) {
        int j = __builtin_ctz(unsigned(mask));
        Split prev = list[mask & (mask-1)];
        int w = f.atoms[j].power;
        list[mask] = {prev.product*w, std::max(w, prev.maximum), mask, prev.count+1};
    }
    return list;
}

inline bool size_possible(const Split& s, int sum) {
    return sum >= 3 && 2LL*s.maximum+1 <= sum &&
           Wide(27)*s.product*s.product <= Wide(sum)*sum*sum;
}

inline int inverse(int a, int b) {
    int64_t r = a, q = b, u = 1, v = 0;
    while (q) {
        int64_t k = r/q, nr = r-k*q, nu = u-k*v;
        r = q; q = nr; u = v; v = nu;
    }
    if (r != 1) std::abort();
    return int((u % b + b) % b);
}
