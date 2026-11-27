// Exhaust finite carry graphs at fixed k, allowing every positive odd scalar Q.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

using Word = std::uint64_t;
using State = std::array<Word, 3>;
using Wide = __uint128_t;

struct Hash {
    std::size_t operator()(const State &state) const {
        Word h = state[0] * 0x9e3779b97f4a7c15ULL;
        h ^= state[1] + (h << 6) + (h >> 2);
        return h ^ (state[2] + (h << 6) + (h >> 2));
    }
};

enum class Result { rejected, accepted, unfinished };

struct Search {
    std::unordered_set<State, Hash> seen;
    std::vector<State> queue;

    void insert(const State &state) {
        if (seen.insert(state).second) queue.push_back(state);
    }
};

// u holds the pending k+1 bits of Q. f holds the upper quotient of c*Q-1.
// The third entry stores the previous input bit and carry for Q+2Q.
Result binary(Word P, Word c, Word q, Word prev, Word add, std::size_t cap) {
    const Word B = 2 * P;
    Search search;
    search.insert({q, Word((Wide(c)*q-1)/B), 2*prev+add});
    for (std::size_t i = 0; i < search.queue.size(); ++i) {
        auto [u, f, state] = search.queue[i];
        prev = state / 2;
        add = state % 2;
        if (!(f & ~u) && !(prev && (prev+add)%2)) return Result::accepted;
        if (search.queue.size() >= cap) return Result::unfinished;
        for (Word b : {0, 1}) {
            Word sum = b+prev+add;
            if (prev && sum%2) continue;
            Word next = f+b*c;
            if (next%2 > u%2) continue;
            search.insert({u/2+b*P, next/2, 2*b+sum/2});
        }
    }
    return Result::rejected;
}

bool carry_free(State x, Word p) {
    while (x[0] || x[1] || x[2]) {
        if (x[0]%p+x[1]%p+x[2]%p >= p) return false;
        for (Word &a : x) a /= p;
    }
    return true;
}

Result odd_base(Word P, Word c, Word p, std::size_t cap) {
    State slopes = {3*P, c, 2*P-c};
    Search search;
    Word modulus = p==3 ? 9 : 5;
    for (Word q = 0; q < modulus; ++q) {
        if (p==3 && ((5*P*q+1)%3 || (5*P*q+1)%9==0)) continue;
        if (p==5 && P*q%5!=2) continue;
        State upper, lower;
        for (int j = 0; j < 3; ++j) {
            Word part = slopes[j]*q + (j==2);
            upper[j] = part/modulus;
            lower[j] = part%modulus;
        }
        if (carry_free(lower,p)) search.insert(upper);
    }
    for (std::size_t i = 0; i < search.queue.size(); ++i) {
        State s = search.queue[i];
        if (carry_free(s,p)) return Result::accepted;
        if (search.queue.size() >= cap) return Result::unfinished;
        for (Word b = 0; b < p; ++b) {
            State next;
            Word sum = 0;
            for (int j = 0; j < 3; ++j) {
                Word n = s[j]+b*slopes[j];
                sum += n%p;
                next[j] = n/p;
            }
            if (sum<p) search.insert(next);
        }
    }
    return Result::rejected;
}

struct Exponent {
    int k;
    Word P;
    bool witnesses;
    Word tested = 0, failed = 0;
    std::array<Word,3> excluded{};

    void check(Word q, Word prev, Word add) {
        // Newton lifting doubles the inverse precision six times from one bit.
        Word inverse = 1;
        for (int j = 0; j < 6; ++j) inverse *= 2-q*inverse;
        Word c = inverse % (2*P);
        if (c%15) return;
        ++tested;
        std::array<bool,3> accepted{};
        for (std::size_t cap : {10000, 100000, 2200000}) {
            for (int i = 0; i < 3; ++i) {
                if (accepted[i]) continue;
                auto result = i==0 ? binary(P,c,q,prev,add,cap)
                                   : odd_base(P,c,i==1 ? 3 : 5,cap);
                if (result==Result::rejected) {
                    ++excluded[i];
                    if (witnesses) std::cout << "excluded " << k << " " << c
                                             << " " << (i==0 ? 2 : i==1 ? 3 : 5) << '\n';
                    return;
                }
                accepted[i] = result==Result::accepted;
            }
        }
        ++failed;
        std::cout << "unresolved " << k << " " << c << " " << q << std::endl;
    }

    void prefixes(int j, Word q, Word prev, Word add) {
        if (j==k+1) {
            check(q,prev,add);
            return;
        }
        for (Word b : {0, 1}) {
            Word sum = b+prev+add;
            if (prev && sum%2) continue;
            prefixes(j+1,q+(b<<j),b,sum/2);
        }
    }
};

int main(int argc, char **argv) {
    int first = argc>1 ? std::atoi(argv[1]) : 3;
    int last = argc>2 ? std::atoi(argv[2]) : 27;
    bool witnesses = argc==4 && std::string(argv[3])=="--witnesses";
    if (argc>4 || (argc==4 && !witnesses) || first<1 || last<first || last>45) return 2;
    for (int k = first; k <= last; ++k) {
        Exponent run{k, Word(1)<<k, witnesses};
        run.prefixes(1,1,1,0);
        std::cout << k << " " << run.tested;
        for (Word n : run.excluded) std::cout << " " << n;
        std::cout << " " << run.failed << std::endl;
    }
    return 0;
}
