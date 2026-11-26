// Enumerate fixed full prime powers and joint binary words of arbitrary length.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using State = std::array<int,6>;

bool disjoint_low(const State& entries, int mask) {
    for (int i : {0,3}) {
        int a = entries[i], b = entries[i+1], c = entries[i+2];
        if (((a & b) | (a & c) | (b & c)) & mask) return false;
    }
    return true;
}

bool accepts_some_word(const State& slopes, const std::vector<int>& low_digits) {
    // Carries after reading a candidate low byte. Zeros above the final bit
    // never alter acceptance, so words with fewer than eight digits are included.
    std::vector<State> states;
    auto add = [&](const State& state) {
        if (std::find(states.begin(),states.end(),state) == states.end())
            states.push_back(state);
    };
    for (int q : low_digits) {
        State current;
        for (int i = 0; i < 6; ++i) current[i] = slopes[i]*q + (i%3 == 2);
        if (!disjoint_low(current,255)) continue;
        for (int& c : current) c /= 256;
        add(current);
    }
    for (std::size_t k = 0; k < states.size(); ++k) {
        State current = states[k];
        if (disjoint_low(current,-1)) return true;
        for (int bit : {0,1}) {
            State next;
            for (int i = 0; i < 6; ++i) next[i] = current[i]+bit*slopes[i];
            if (!disjoint_low(next,1)) continue;
            for (int& c : next) c /= 2;
            add(next);
        }
    }
    return false;
}

std::vector<int> prime_powers(int bound) {
    std::vector<bool> composite(bound+1,false);
    std::vector<int> powers;
    for (int p = 2; p <= bound; ++p) {
        if (composite[p]) continue;
        for (int k = 2*p; k <= bound; k += p) composite[k] = true;
        for (std::int64_t q = p; q <= bound; q *= p) powers.push_back(static_cast<int>(q));
    }
    std::sort(powers.begin(),powers.end());
    return powers;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int bound = std::stoi(argv[1]);
    // With P <= 500,000, all six slopes times a low byte fit in signed int.
    if (bound < 2 || bound > 500000) return 2;
    std::uint64_t tested7=0, tested9=0, accepted7=0, accepted9=0;
    for (int p : prime_powers(bound)) {
        for (int v : {7,9}) {
            if (std::gcd(p, v == 7 ? 42 : 1584) != 1) continue;
            std::vector<int> low;
            for (int q = 1; q < 256; q += 2) {
                int t = p*q;
                if (t%16 != (v == 7 ? 11 : 3)) continue;
                int a = v == 7 ? 3 : 11, b = v == 7 ? 6 : 4;
                // Only the low byte matters in each cut.
                if ((a*t & ((16-a)*t+1) & 255) ||
                    (b*t & ((v-b)*t+1) & 255)) continue;
                low.push_back(q);
            }
            auto test = [&](int c, int e) {
                if (std::gcd(c,p) != 1) return;
                State slopes = v == 7 ? State{6*p,c,p-c,3*p,e,13*p-e}
                                      : State{4*p,c,5*p-c,11*p,c,5*p-c};
                (v == 7 ? tested7 : tested9)++;
                if (accepts_some_word(slopes,low)) {
                    (v == 7 ? accepted7 : accepted9)++;
                    std::cout << p << ' ' << v << ' ' << c << ' ' << e << '\n';
                }
            };
            if (v == 7) {
                for (int c = 336; c < p; c += 336) test(c,c);
                for (int c = (84*p)%336; c < p; c += 336) test(c,c+12*p);
            } else {
                for (int c = 1584; c < 5*p; c += 1584) test(c,c);
            }
        }
    }
    std::cerr << "Fixed power bound " << bound << ", (7,16) shapes " << tested7
              << ", nonempty " << accepted7 << ", (9,16) shapes " << tested9
              << ", nonempty " << accepted9 << '\n';
}
