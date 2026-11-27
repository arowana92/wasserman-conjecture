// Enumerate binary scales for the fixed marker 2^31-1, at every outer exponent.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

using Integer = std::int64_t;
constexpr Integer marker = (Integer(1) << 31) - 1;
using Key = std::tuple<int, int, int>;

struct Shape {
    int a, b;
    std::vector<bool> scale_residues;
};

Integer inverse(Integer a, Integer m) {
    Integer b = m, x = 1, y = 0;
    while (b) {
        Integer q = a / b;
        Integer next_r = a - q * b, next_x = x - q * y;
        a = b; b = next_r; x = y; y = next_x;
    }
    return (x % m + m) % m;
}

std::vector<int> odd_powers(int d) {
    std::vector<int> powers;
    for (int p = 3; p <= d; p += 2) {
        if (d % p) continue;
        int w = 1;
        while (d % p == 0) { w *= p; d /= p; }
        powers.push_back(w);
    }
    return powers;
}

std::vector<bool> allowed_residues(int d, int i, int j, int a, int b) {
    auto powers = odd_powers(d);
    std::vector<bool> good(d);
    for (int u = 1; u < d; ++u) {
        if (std::gcd(u, d) != 1) continue;
        Integer z = 1;
        do {
            bool pass = true;
            for (int w : powers) {
                Integer x = marker * a % w, y = marker * b % w;
                bool second = (!x && !y) || ((!x || !y) && (x + y - u) % w == 0);
                Integer f = i * u * z % w, g = j * u * z % w;
                bool first = (!f && !g) || ((!f || !g) && (f + g) % w == 1);
                if (!(first || second)) { pass = false; break; }
            }
            if (pass) { good[u] = true; break; }
            z = 2 * z % d;
        } while (z != 1);
    }
    return good;
}

void scan(int d, int i, int j, const std::vector<Shape>& shapes) {
    Integer total = d * marker;
    Integer high = (total - 1) / (i + j), low = high;
    for (const auto& shape : shapes)
        low = std::min(low, marker * (shape.a + shape.b) + 1);

    auto visit = [&](auto&& self, Integer r, Integer step) -> void {
        Integer x = i * r, y = j * r;
        Integer mask = step - 1;
        if (((x & y) | ((x | y) & ~total)) & mask) return;
        Integer u = r;
        if (u < low) u += ((low - u + step - 1) / step) * step;
        if (u > high) return;
        if (step <= high - low) {
            self(self, r, 2 * step);
            self(self, r + step, 2 * step);
            return;
        }
        x = i * u; y = j * u;
        if ((x & y) || ((x | y) & ~total) || u % marker == 0) return;
        Integer v = 0;
        for (const auto& shape : shapes) {
            if (u <= marker * (shape.a + shape.b) || !shape.scale_residues[u % d])
                continue;
            if (!v) v = inverse(u % marker, marker);
            Integer z = v * shape.a % marker + v * shape.b % marker;
            std::uint32_t classes = 0;
            for (int e = 0; e < 31; ++e)
                if (z <= (Integer(d) << e) % marker) classes |= 1u << e;
            if (classes)
                std::cout << d << ' ' << u << ' ' << i << ' ' << j << ' '
                          << shape.a << ' ' << shape.b << ' ' << classes << '\n';
        }
    };
    visit(visit, 1, 2);
}

int main() {
    std::map<Key, std::vector<Shape>> groups;
    int d, i, j, a, b;
    while (std::cin >> d >> i >> j >> a >> b) {
        if (d < 3 || d > 511 || !(d & 1) || i < 1 || i >= j || a < 1 || a > b ||
            (i + j) * (a + b) >= d) return 1;
        auto residues = allowed_residues(d, i, j, a, b);
        if (std::count(residues.begin(), residues.end(), true))
            groups[{d, i, j}].push_back({a, b, std::move(residues)});
    }
    for (const auto& entry : groups) {
        const auto [d, i, j] = entry.first;
        scan(d, i, j, entry.second);
    }
}
