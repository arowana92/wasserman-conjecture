// Find a row avoiding each four-element subset of the supplied endpoint primes.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Integer = std::uint64_t;
using Wide = __uint128_t;
using Layer = std::pair<Integer, Integer>;

Integer inverse(Integer value, Integer prime) {
    Integer result = 1;
    for (Integer exponent = prime - 2; exponent; exponent /= 2) {
        if (exponent % 2) result = Wide(result) * value % prime;
        value = Wide(value) * value % prime;
    }
    return result;
}

std::pair<Integer, Integer> find_row(Integer n, const std::vector<Integer>& primes) {
    Wide period = 1;
    std::vector<Integer> residues{0};
    std::vector<Layer> layers;
    for (Integer p : primes) {
        std::vector<Integer> next;
        Integer inv = inverse(period % p, p);
        for (Integer r : residues) {
            for (Integer v = 0; v <= n % p; ++v) {
                Integer step = Wide((p + v - r % p) % p) * inv % p;
                Wide x = r + period * step;
                if (x < n) next.push_back(static_cast<Integer>(x));
            }
        }
        residues.swap(next);
        period *= p;
        for (Integer w = p; w <= n / p;) {
            w *= p;
            layers.emplace_back(w, n % w);
        }
    }

    // If a mod w > n mod w, every intervening cut through the next
    // multiple of w minus one also fails this factorial layer.
    std::vector<Integer> cuts;
    for (Integer residue : residues) {
        Wide current = residue ? residue : period;
        while (current <= n / 2) {
            Integer a = static_cast<Integer>(current);
            Integer skip = 0;
            for (auto [w, limit] : layers) {
                Integer r = a % w;
                if (r > limit) skip = std::max(skip, w - r);
            }
            if (!skip) cuts.push_back(a);
            current += period * (skip ? 1 + (skip - 1) / period : 1);
        }
    }

    std::sort(cuts.begin(), cuts.end());
    for (Integer p : primes) layers.emplace_back(p, n % p);
    for (std::size_t i = 0; i < cuts.size(); ++i) {
        Integer a = cuts[i];
        for (std::size_t j = i; j < cuts.size(); ++j) {
            Integer b = cuts[j];
            if (b > (n - a) / 2) break;
            bool free = std::all_of(layers.begin(), layers.end(), [=](Layer l) {
                auto [w, limit] = l;
                return a % w + b % w <= limit;
            });
            if (free) return {a, b};
        }
    }
    return {0, 0};
}

void check_input(Integer n, const std::vector<Integer>& primes) {
    if (n < 3) throw std::invalid_argument("weight must be at least three");
    for (auto it = primes.begin(); it != primes.end(); ++it) {
        Integer p = *it;
        if (p < 2 || p > 1000000 || n % p > 1 ||
            std::find(primes.begin(), it, p) != it)
            throw std::invalid_argument("expected distinct endpoint primes at most one million");
        for (Integer j = 2; j <= p / j; ++j)
            if (p % j == 0) throw std::invalid_argument("composite base");
    }
}

void print_row(Integer n, const std::vector<Integer>& group) {
    auto [a, b] = find_row(n, group);
    for (Integer p : group) std::cout << p << ' ';
    std::cout << a << ' ' << b << '\n';
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--groups") {
        // Each line contains a weight followed by one group of endpoint primes.
        std::string line;
        while (std::getline(std::cin, line)) {
            std::istringstream input(line);
            Integer n, p;
            if (!(input >> n)) throw std::invalid_argument("missing weight");
            std::vector<Integer> primes;
            while (input >> p) primes.push_back(p);
            if (!input.eof() || primes.size() < 2 || primes.size() > 5)
                throw std::invalid_argument("supply between two and five primes");
            check_input(n, primes);
            print_row(n, primes);
        }
        return 0;
    }
    if (argc < 6) throw std::invalid_argument("supply the weight and at least four primes");
    Integer n = std::stoull(argv[1]);
    std::vector<Integer> primes;
    for (int i = 2; i < argc; ++i) primes.push_back(std::stoull(argv[i]));
    check_input(n, primes);
    for (std::size_t i = 0; i < primes.size(); ++i)
        for (std::size_t j = i + 1; j < primes.size(); ++j)
            for (std::size_t k = j + 1; k < primes.size(); ++k)
                for (std::size_t l = k + 1; l < primes.size(); ++l) {
                    std::vector<Integer> group{primes[i], primes[j], primes[k], primes[l]};
                    print_row(n, group);
                }
}
