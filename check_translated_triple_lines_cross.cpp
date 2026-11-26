// Independent upper bounds for the two rational gcd degrees by reduction modulo primes.
// Compile with c++ -O3 -std=c++17 check_translated_triple_lines_cross.cpp.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using Vec = std::array<int, 3>;
using Grid = std::array<int, 6>;
using Poly = std::vector<int>;

Grid positions(Vec v) {
    return {0, v[0], v[1], v[2], v[2]+v[0], v[2]+v[1]};
}

bool triangle(int a, int b, int c, int p) {
    Vec v = {a, b, c};
    std::sort(v.begin(), v.end());
    int gap = v[2]-v[0];
    return v[0] == v[1] && (gap == 2 || (p == 2 && gap == 4));
}

bool small_offset(Vec u, int type, int p) {
    auto [x,y,z] = u;
    if (type == 0)
        return (x == y+z && (triangle(x+z,y,z,p) || triangle(0,x,y,p)))
            || (z == 0 && triangle(0,x,y+z,p));
    if (type == 1)
        return (x == z && (triangle(x+z,y,y+z,p) || triangle(0,x,y,p)))
            || (y == z && triangle(0,x,z,p));
    return (x == 0 || y+z == 0) && triangle(x,y,y+z,p);
}

bool feasible(Vec v, Vec e, int type, int p) {
    Grid a = positions(v), b = positions(e);
    auto same = [&](int i, int j) { return a[i] == a[j] && b[i] == b[j]; };
    if (same(0,1) || same(0,2) || same(1,2) || same(0,3)) return false;
    if ((type == 1 && same(2,3)) || (type == 2 && same(5,0))) return false;
    for (int i = 0; i < 6 && p == 2; ++i)
        for (int j = 0; j < i; ++j)
            if (same(i,j)) return false;
    return true;
}

int power(int a, int b, int q) {
    int t = 1;
    while (b) {
        if (b % 2) t = static_cast<long long>(t)*a % q;
        a = static_cast<long long>(a)*a % q;
        b /= 2;
    }
    return t;
}

void trim(Poly &f) {
    while (!f.empty() && f.back() == 0) f.pop_back();
}

int gcd_degree(Poly f, Poly g, int q) {
    while (!g.empty()) {
        int inv = power(g.back(), q-2, q);
        while (f.size() >= g.size()) {
            int coeff = static_cast<long long>(q-f.back())*inv % q;
            int offset = static_cast<int>(f.size()-g.size());
            for (int i = 0; i < static_cast<int>(g.size()); ++i)
                f[offset+i] = (f[offset+i]+static_cast<long long>(coeff)*g[i]) % q;
            trim(f);
        }
        f.swap(g);
    }
    return static_cast<int>(f.size())-1;
}

bool certifies(Vec v, Vec e, int type, int p, int q) {
    Grid a = positions(v), b = positions(e);
    int shift = *std::min_element(b.begin(), b.end());
    for (int &c : b) c = power(p, c-shift, q);
    auto sum = [&](const std::vector<int> &indices) {
        int low = a[indices[0]], high = low;
        for (int i : indices) { low = std::min(low,a[i]); high = std::max(high,a[i]); }
        Poly f(high-low+1);
        for (int i : indices) f[a[i]-low] = (f[a[i]-low]+b[i]) % q;
        return f;
    };
    const std::array<std::vector<int>,3> pairs = {{{1,5},{1,3},{3,4}}};
    const std::array<std::vector<int>,3> triples = {{{2,3,4},{2,4,5},{1,2,5}}};
    Poly f = sum(pairs[type]), g = sum(triples[type]), m = sum({0,1,2,3,4,5});
    if (!f.back() || !g.back() || !m.back()) return false;
    int h = static_cast<int>(m.size())-1;
    Poly fg(f.size()+g.size()-1);
    for (size_t i = 0; i < f.size(); ++i)
        for (size_t j = 0; j < g.size(); ++j)
            fg[i+j] = (fg[i+j]+static_cast<long long>(f[i])*g[j]) % q;
    return 3*(2*h-3*gcd_degree(f,g,q)-2*gcd_degree(m,fg,q)) >= h;
}

int main() {
    long long checked = 0, second_modulus = 0;
    for (int p : {2,3}) for (int type = 0; type < 3; ++type) {
        std::vector<Vec> offsets;
        for (int x = -8; x <= 8; ++x) for (int y = -8; y <= 8; ++y)
            for (int z = -8; z <= 8; ++z)
                if (small_offset({x,y,z},type,p)) offsets.push_back({x,y,z});
        const int count[2][3] = {{13,14,10},{7,7,5}};
        assert(static_cast<int>(offsets.size()) == count[p-2][type]);
        int limit = type == 0 ? 36 : 21;
        for (int x = 0; x <= limit; ++x) for (int y = -limit; y <= limit; ++y)
            for (int z = -limit; z <= limit; ++z) {
                if (x == 0 && (y < 0 || (y == 0 && z < 0))) continue;
                if (std::gcd(x,std::gcd(y,z)) != 1) continue;
                if (std::max({0,x,y})-std::min({0,x,y})+std::abs(z) > limit) continue;
                Vec v = {x,y,z};
                for (Vec t : offsets) for (int s = 0; s < 4; ++s) {
                    Vec e;
                    bool integral = true, parallel = true;
                    for (int i = 0; i < 3; ++i) {
                        int n = t[i]+s*v[i];
                        if (n % 4) integral = false;
                        e[i] = n/4;
                        int j = (i+1)%3;
                        if (t[i]*v[j] != t[j]*v[i]) parallel = false;
                    }
                    if (!integral || parallel || !feasible(v,e,type,p)) continue;
                    ++checked;
                    if (certifies(v,e,type,p,1009)) continue;
                    ++second_modulus;
                    if (certifies(v,e,type,p,10007) || certifies(v,e,type,p,65537)) continue;
                    std::cerr << "Uncertified line at p=" << p << " type=" << type
                              << " v=" << x << ',' << y << ',' << z << '\n';
                    return 1;
                }
            }
    }
    std::cout << checked << " incidences certified, " << second_modulus
              << " needed a second modulus\n";
}
