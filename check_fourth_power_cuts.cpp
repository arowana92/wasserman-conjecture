// Enumerate positive odd s < 3^M, 3 coprime to s, for which 16s^4
// and 19s^4 add without a ternary carry.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

using Int = std::uint64_t;
using Wide = unsigned __int128;
constexpr int chunk = 6561; // 3^8
constexpr Int limb = Int(chunk) * chunk * chunk; // 3^24

struct TritMask { unsigned nonzero, two; };
std::array<TritMask, chunk> masks;

void initialize_masks() {
    for (int v = 1; v < chunk; ++v) {
        auto upper = masks[v / 3];
        masks[v] = {(upper.nonzero << 1) | unsigned(v % 3 != 0),
                    (upper.two << 1) | unsigned(v % 3 == 2)};
    }
}

bool chunk_clear(int a, int b) {
    auto x = masks[a], y = masks[b];
    return !(x.nonzero & y.two) && !(y.nonzero & x.two);
}

template <class T>
bool clear(T a, T b) {
    while (a && b) {
        if (!chunk_clear(int(a % chunk), int(b % chunk))) return false;
        a /= chunk;
        b /= chunk;
    }
    return true;
}

Int pow3(int n) {
    Int v = 1;
    for (int k = 0; k < n; ++k) v *= 3;
    return v;
}

Int fourth_mod(Int s, Int m) {
    Int square = Wide(s) * s % m;
    return Wide(square) * square % m;
}

std::vector<Int> low_classes(int digits) {
    std::vector<Int> current = {1, 2};
    Int place = 3;
    for (int k = 1; k < digits; ++k) {
        std::vector<Int> next;
        next.reserve(2 * current.size());
        for (Int d = 0; d < 3; ++d) {
            for (Int v : current) {
                Int s = d * place + v, w = fourth_mod(s, 3 * place);
                if ((16 * w / place) % 3 + (19 * w / place) % 3 < 3)
                    next.push_back(s);
            }
        }
        assert(next.size() == 2 * current.size());
        current.swap(next);
        place *= 3;
    }
    return current;
}

// h <= 3^19: both fourth powers times 19 fit in 128 bits.
Wide fourth(Int h) {
    Wide square = Wide(h) * h;
    return square * square;
}

bool possible_leading_prefix(Int h, Wide scale) {
    Wide left = fourth(h), upper = fourth(h + 1);
    for (; scale < 19 * (upper - left); scale *= 3) {}
    Wide x0 = 16 * left / scale, x1 = (16 * upper - 1) / scale;
    Wide y0 = 19 * left / scale, y1 = (19 * upper - 1) / scale;
    for (Int x = Int(x0); x <= x1; ++x)
        for (Int y = Int(y0); y <= y1; ++y)
            if (clear(x, y)) return true;
    return false;
}

std::vector<Int> high_prefixes(int digits) {
    std::vector<Int> current = {0};
    Wide nonzero_scale = 1;
    Int leading_place = 1;
    for (int k = 1; k <= digits; ++k) {
        std::vector<Int> next;
        for (Int v : current) {
            for (Int d = 0; d < 3; ++d) {
                Int h = 3 * v + d;
                // For smaller prefixes (including leading zeroes) use the
                // actual size of the interval to pick the digit boundary.
                Wide scale = nonzero_scale;
                for (Int place = leading_place; place > h && place > 1; place /= 3)
                    scale /= 27;
                if (possible_leading_prefix(h, scale)) next.push_back(h);
            }
        }
        current.swap(next);
        nonzero_scale *= 27;
        leading_place *= 3;
    }
    return current;
}

// Use base 3^24 for the exact fourth power. For M <= 38 the square of s
// fits in 128 bits and needs at most four limbs. The convolution fits too.
bool full_test(Int s) {
    Wide square = Wide(s) * s;
    std::array<Int, 4> digits{};
    int n = 0;
    while (square) {
        digits[n++] = square % limb;
        square /= limb;
    }
    std::array<Wide, 9> product{};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            product[i+j] += Wide(digits[i]) * digits[j];
    Int carry_x = 0, carry_y = 0;
    for (int k = 0; k < 2*n; ++k) {
        product[k+1] += product[k] / limb;
        Int digit = product[k] % limb;
        Int x = 16 * digit + carry_x;
        Int y = 19 * digit + carry_y;
        if (!clear(x % limb, y % limb)) return false;
        carry_x = x / limb;
        carry_y = y / limb;
    }
    assert(product[2*n] == 0);
    return clear(carry_x, carry_y);
}

int main(int argc, char** argv) {
    int digits = argc > 1 ? std::atoi(argv[1]) : 32;
    int head = argc > 2 ? std::atoi(argv[2]) : std::min(19, (digits*2)/3);
    int tail = argc > 3 ? std::atoi(argv[3]) : std::min(23, (digits*3)/4);
    if (digits < 3 || digits > 38 || head < 1 || head > 19 || head > digits ||
        tail < 1 || tail > 27 || tail > digits) throw std::invalid_argument("digit bounds");
    initialize_masks();
    auto low = low_classes(tail);
    auto high = high_prefixes(head);
    Int interval = pow3(digits-head), modulus = pow3(tail);
    Int joined = 0, tested = 0, found = 0;
    const Int middle = pow3(36);
    for (Int h : high) {
        Int left = h*interval, end = left+interval;
        for (Int base = (left/modulus)*modulus; base < end; base += modulus) {
            auto it = std::lower_bound(low.begin(), low.end(), left > base ? left-base : 0);
            while (it != low.end() && *it < end-base) {
                Int s = base + *it++;
                if (!(s & 1)) continue;
                ++joined;
                Int w = fourth_mod(s, middle);
                if (!clear((16*w) % middle, (19*w) % middle)) continue;
                ++tested;
                if (full_test(s)) {
                    ++found;
                    std::cout << s << '\n';
                }
            }
        }
    }
    std::cerr << "high " << high.size() << " low " << low.size()
              << " joined " << joined << " tested " << tested << " clear " << found << '\n';
}
