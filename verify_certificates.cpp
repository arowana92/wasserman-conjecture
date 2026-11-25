// Compile with c++ -std=c++17 -O3 verify_certificates.cpp -o verify_certificates
// Usage: ./verify_certificates FIRST LAST < certificate-file
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Mask = unsigned;

bool is_prime(int p) {
    if (p < 2) return false;
    for (int d = 2; d <= p / d; ++d)
        if (p % d == 0) return false;
    return true;
}

using Valuations = std::vector<std::vector<int>>;

// Legendre's formula, summed by the recurrence v_p(i!) = floor(i/p) + v_p(floor(i/p)!).
Valuations valuations(int n, const std::vector<int>& primes) {
    Valuations v(primes.size(), std::vector<int>(n + 1));
    for (size_t j = 0; j < primes.size(); ++j)
        for (int i = 1; i <= n; ++i) {
            const int quotient = i / primes[j];
            v[j][i] = quotient + v[j][quotient];
        }
    return v;
}

std::map<Mask, std::uint64_t> histogram(int n, const std::vector<int>& primes) {
    auto v = valuations(n, primes);
    std::vector<std::uint64_t> counts(1u << primes.size());
    std::vector<Mask> masks(n + 1);
    for (int a = 1; a <= n / 3; ++a) {
        const int last = (n - a) / 2;
        std::fill(masks.begin() + a, masks.begin() + last + 1, 0);
        for (size_t j = 0; j < primes.size(); ++j) {
            const auto& col = v[j];
            const int remainder = col[n] - col[a];
            for (int b = a; b <= last; ++b)
                if (remainder > col[b] + col[n-a-b]) masks[b] |= 1u << j;
        }
        for (int b = a; b <= last; ++b) ++counts[masks[b]];
    }
    std::map<Mask, std::uint64_t> present;
    for (Mask mask = 0; mask < counts.size(); ++mask)
        if (counts[mask]) present[mask] = counts[mask];
    return present;
}

int main(int argc, char** argv) {
    if (argc != 3) throw std::runtime_error("Specify the first and last weight");
    int n = std::stoi(argv[1]);
    const int first = n, last = std::stoi(argv[2]);
    if (first < 3 || last < first || last > 1000000)
        throw std::runtime_error("Require 3 <= FIRST <= LAST <= 1000000");
    std::string line;
    std::uint64_t total = 0;
    while (std::getline(std::cin, line)) {
        std::istringstream input(line);
        int listed, p;
        if (!(input >> listed) || listed != n || listed > last)
            throw std::runtime_error("Missing, extra or out-of-order weight");
        std::vector<int> primes;
        while (input >> p) {
            if (p > n || !is_prime(p) || (!primes.empty() && primes.back() >= p))
                throw std::runtime_error("Expected increasing primes at " + std::to_string(n));
            primes.push_back(p);
        }
        if (!input.eof() || primes.empty() || primes.size() > 12)
            throw std::runtime_error("Malformed certificate at " + std::to_string(n));
        auto seen = histogram(n, primes);
        for (auto [x, count] : seen) {
            total += count;
            for (auto [y, cy] : seen)
                for (auto [z, cz] : seen)
                    if (!(x & y & z)) throw std::runtime_error("Failed certificate at " + std::to_string(n));
            if (first == last) {
                std::cout << "mask " << x << " count " << count << "\n";
            }
        }
        ++n;
    }
    if (n != last + 1) throw std::runtime_error("Incomplete interval");
    std::cout << "Verified " << first << ".." << last << ", " << total
              << " unordered positive triples.\n";
}
