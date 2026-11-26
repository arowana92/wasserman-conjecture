// Enumerate quadratic-position allocations left by the rational inequalities.
// Input comes from check_quadratic_digit_positions.py, which also checks the output.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

using DegreeTable = std::map<std::tuple<int,int,int,int>, std::map<int,int>>;
constexpr int height_cap = 68333;

struct Direction {
    int q, r, v;
    std::map<int,int> exceptional_degree;
};

std::vector<Direction> directions(int i, int j, const DegreeTable& degrees) {
    std::vector<Direction> result;
    for (int v = 1; v <= 100; ++v) {
        for (int r = 0; (i+j)*r <= v; ++r) {
            int q = v-(i+j)*r;
            if (std::gcd(q,r) != 1) continue;
            auto found = degrees.find({i,j,v,r});
            result.push_back({q,r,v,found == degrees.end() ? std::map<int,int>{} : found->second});
        }
    }
    return result;
}

// Skip previously certified values of t. Each newly covered point gets
// its own evaluation of the original inequality using integer arithmetic.
int next_uncovered(std::vector<int>& next, int t) {
    while (next[t] != t) {
        next[t] = next[next[t]];
        t = next[t];
    }
    return t;
}

void mark_interval(std::vector<int>& next, int c, int i, int j, int error, const Direction& dir) {
    int center = i+j, gap = j-i;
    int t_weight = 10-6*(center == 5), c_weight = 50-6*center*(center == 5);
    // Intersect the two linear half-lines from
    // 27 gap v + 3 error |c q - t r| <= t_weight t + c_weight c.
    int lower = 1, upper = static_cast<int>(next.size())-2;
    int rhs = 27*gap*dir.v + (3*error*dir.q-c_weight)*c;
    int coeff = 3*error*dir.r+t_weight;
    if (rhs > 0) lower = std::max(lower, 1+(rhs-1)/coeff);
    coeff = 3*error*dir.r-t_weight;
    rhs = (3*error*dir.q+c_weight)*c-27*gap*dir.v;
    if (coeff > 0) {
        if (rhs < 0) return;
        upper = std::min(upper, rhs/coeff);
    } else if (rhs < 0) {
        if (!coeff) return;
        lower = std::max(lower, 1+(-rhs-1)/(-coeff));
    }
    if (lower > upper) return;
    for (int t = next_uncovered(next,lower); t <= upper; t = next_uncovered(next,t+1)) {
        int delta = c*dir.q-t*dir.r;
        if (!delta) continue;
        int degree = 2*dir.v*(center == 5);
        auto found = dir.exceptional_degree.find(delta);
        if (found != dir.exceptional_degree.end()) degree = found->second;
        int L = t+center*c, H = 5*t+25*c;
        long long bound = 3LL*degree*L + 3LL*(9*gap*dir.v+error*std::abs(delta))*dir.v;
        if (bound <= 2LL*H*dir.v) next[t] = next_uncovered(next,t+1);
    }
}

void check_pair(int i, int j, const DegreeTable& degrees) {
    int center = i+j;
    std::vector<int> heights;
    for (int k = 0; k < 6; ++k) {
        if (k != i && k != j) heights.push_back(k*(k-center));
    }
    auto [low, high] = std::minmax_element(heights.begin(), heights.end());
    int error = 5*(j-i)*(*high-*low);
    auto slopes = directions(i,j,degrees);
    for (int c = 1; 25*c < height_cap; ++c) {
        int max_t = (height_cap-25*c)/5;
        std::vector<int> next(max_t+2);
        std::iota(next.begin(), next.end(), 0);
        for (const auto& dir : slopes) mark_interval(next, c, i, j, error, dir);
        for (int t = 1; t <= max_t; ++t) {
            if (next[t] == t && std::gcd(t,c) < 5) {
                std::cout << i << ' ' << j << ' ' << t << ' ' << c << '\n';
            }
        }
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    DegreeTable degrees;
    int i, j, v, r, delta, degree;
    while (std::cin >> i >> j >> v >> r >> delta >> degree) {
        degrees[{i,j,v,r}][delta] = degree;
    }
    for (int a = 0; a < 6; ++a) {
        for (int b = a+1; b < 6; ++b) check_pair(a,b,degrees);
    }
}
