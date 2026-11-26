// Enumerate the necessary row when N-1 has three prime-power factors.
// Usage: check_three_factors [largest prime power [--cut-pairs]], default 20000.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using Integer = std::int64_t;
using Row = std::array<Integer, 3>;

struct PrimePower {
    int prime;
    int power;
};

std::vector<PrimePower> prime_powers(int bound) {
    std::vector<PrimePower> result;
    std::vector<bool> composite(bound + 1);
    for (int p = 2; p <= bound; ++p) {
        if (composite[p]) continue;
        for (Integer j = Integer(p)*p; j <= bound; j += p) composite[j] = true;
        for (Integer w = p; w <= bound; w *= p) result.push_back({p, int(w)});
    }
    std::sort(result.begin(), result.end(), [](PrimePower a, PrimePower b) {
        return a.power < b.power;
    });
    return result;
}

int inverse(int value, int modulus) {
    Integer a = value, b = modulus, x = 1, y = 0;
    while (b != 0) {
        Integer q = a / b;
        a -= q*b;
        std::swap(a, b);
        x -= q*y;
        std::swap(x, y);
    }
    if (a != 1) throw std::logic_error("moduli must be coprime");
    return int((x % modulus + modulus) % modulus);
}

bool avoids(Row row, int prime) {
    auto [a, b, c] = row;
    if (prime == 2) return !(a & b) && !(b & c) && !(a & c);
    while (a || b || c) {
        if (a % prime + b % prime + c % prime >= prime) return false;
        a /= prime;
        b /= prime;
        c /= prime;
    }
    return true;
}

template <typename Predicate>
Integer avoided_part_of_weight(Integer n, int common_factor, Predicate avoids_base) {
    Integer product = 1;
    auto include = [&](int p) {
        if (avoids_base(p)) {
            while (n % p == 0) {
                product *= p;
                n /= p;
            }
        }
    };
    for (int p = 2; p <= common_factor/p; ++p) {
        if (common_factor % p != 0) continue;
        include(p);
        do { common_factor /= p; } while (common_factor % p == 0);
    }
    if (common_factor > 1) include(common_factor);
    return product;
}

bool cut_avoids(Integer n, Integer a, int p) {
    if (p == 2) return !(a & (n-a));
    while (a) {
        if (a % p > n % p) return false;
        a /= p;
        n /= p;
    }
    return true;
}

struct PairCounts {
    Integer two_cuts = 0;
    Integer candidates = 0;
};

void check_pairs(Integer n, Row residues, const std::array<PrimePower, 3>& atoms,
                 PairCounts& counts) {
    if (residues[2] < 0) residues[2] += n-1;
    std::array<bool, 3> good;
    for (int i = 0; i < 3; ++i) {
        good[i] = std::all_of(atoms.begin(), atoms.end(), [&](PrimePower atom) {
            return cut_avoids(n, residues[i], atom.prime);
        });
    }
    if (std::count(good.begin(), good.end(), true) < 2) return;
    ++counts.two_cuts;
    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            if (!good[i] || !good[j]) continue;
            Integer x = residues[i], y = residues[j];
            // gcd(n,x) divides the inverse multiplier at i and fits in int.
            int g = static_cast<int>(std::gcd(std::gcd(n, x), y));
            if (g < 3) continue;
            Integer u = avoided_part_of_weight(n, g, [&](int p) {
                return cut_avoids(n, x, p) && cut_avoids(n, y, p);
            });
            if (u < 3) continue;
            ++counts.candidates;
            std::cout << n << ' ' << u << ' ' << x << ' ' << y;
            for (auto atom : atoms) std::cout << ' ' << atom.power;
            std::cout << '\n';
        }
    }
}

void check(int bound, bool cut_pairs) {
    auto atoms = prime_powers(bound);
    int size = static_cast<int>(atoms.size());
    // inv[i][j] is the inverse of power j modulo power i, or 0 for shared bases.
    std::vector<std::vector<int>> inv(size, std::vector<int>(size));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (atoms[i].prime != atoms[j].prime)
                inv[i][j] = inverse(atoms[j].power, atoms[i].power);
        }
    }

    PairCounts pairs;
    Integer weights = 0, correct_sum = 0, common_factor = 0, avoids_shift = 0, eligible = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            if (!inv[i][j]) continue;
            Integer a = atoms[i].power, b = atoms[j].power;
            for (int k = j + 1; k < size; ++k) {
                if (!inv[i][k] || !inv[j][k]) continue;
                ++weights;
                Integer c = atoms[k].power, n = a*b*c + 1;
                Integer x = Integer(inv[i][j])*inv[i][k] % a;
                Integer y = Integer(inv[j][i])*inv[j][k] % b;
                Row row = {b*c*x, a*c*y, n - b*c*x - a*c*y};
                if (cut_pairs) {
                    check_pairs(n, row, {atoms[i], atoms[j], atoms[k]}, pairs);
                    continue;
                }
                // The third value has the prescribed residues. Its positivity
                // puts it in the required interval (0,N-1).
                if (row[2] <= 0) continue;
                ++correct_sum;
                // A common divisor of all parts divides n, hence is coprime to abc.
                int g = int(std::gcd(std::gcd(x, y), row[2]));
                if (g < 3) continue;
                ++common_factor;
                if (!avoids(row, atoms[i].prime) || !avoids(row, atoms[j].prime)
                    || !avoids(row, atoms[k].prime)) continue;
                ++avoids_shift;
                Integer u = avoided_part_of_weight(n, g, [&](int p) { return avoids(row, p); });
                if (u >= 3) {
                    ++eligible;
                    std::cout << "Eligible N " << n << " row " << row[0] << ' '
                              << row[1] << ' ' << row[2] << " U " << u << '\n';
                }
            }
        }
    }
    if (cut_pairs) {
        std::cerr << "weights " << weights << "\nweights with two predecessor-admissible cuts "
                  << pairs.two_cuts << "\ncandidate pairs " << pairs.candidates << '\n';
        return;
    }
    std::cout << "prime-power bound " << bound << "\nweights " << weights
              << "\ncorrect sum " << correct_sum << "\ngcd >= 3 " << common_factor
              << "\navoids shifted primes " << avoids_shift
              << "\neligible " << eligible << '\n';
}

int main(int argc, char** argv) {
    // The bound keeps products of three powers inside a signed 64-bit integer.
    int bound = argc == 1 ? 20000 : std::stoi(argv[1]);
    bool cut_pairs = argc == 3 && std::string(argv[2]) == "--cut-pairs";
    if (argc > 3 || (argc == 3 && !cut_pairs) || bound < 2 || bound > 2000000) {
        std::cerr << "Supply one prime-power bound from 2 to 2000000, optionally followed by --cut-pairs.\n";
        return 2;
    }
    check(bound, cut_pairs);
}
