// Verify the streamed low-digit certificates with a separate prefix generator
// and an automaton that follows the total and two parts by ordinary division.
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using Word = std::uint64_t;
using Wide = __uint128_t;
using State = std::array<Word, 5>;

struct Hash {
    std::size_t operator()(const State& s) const {
        std::size_t value = 0;
        for (Word x : s) value ^= x + 0x9e3779b97f4a7c15ULL + (value<<6) + (value>>2);
        return value;
    }
};

void require(bool ok) {
    if (!ok) throw std::runtime_error("inconsistent certificate");
}

bool fits(Word n, Word a, Word b, Word p) {
    if (p == 2) return !(a & b) && !((a | b) & ~n);
    while (n || a || b) {
        if (a%p + b%p > n%p) return false;
        n /= p;
        a /= p;
        b /= p;
    }
    return true;
}

bool excluded(Word P, Word c, Word q, Word p) {
    Word m = p == 2 ? 2*P : (p == 3 ? 9 : 5);
    std::vector<State> queue;
    std::unordered_set<State, Hash> seen;
    auto insert = [&](const State& s) {
        if (seen.insert(s).second) queue.push_back(s);
    };
    for (Word v = p==2 ? q : 0; v<m; v += p==2 ? m : 1) {
        if (p==3 && ((5*P*v+1)%3 || (5*P*v+1)%9 == 0)) continue;
        if (p==5 && P*v%5 != 2) continue;
        Wide n = Wide(5)*P*v+1, a = Wide(3)*P*v, b = Wide(c)*v;
        if (fits(n%m,a%m,b%m,p))
            insert({Word(n/m),Word(a/m),Word(b/m),p==2 ? 3*v/m : 0,p==2 ? 2*v/m : 0});
    }
    for (std::size_t i=0; i<queue.size(); ++i) {
        auto [n,a,b,f,g] = queue[i];
        if (fits(n,a,b,p) && !(f & g)) return false;
        for (Word bit=0; bit<p; ++bit) {
            Word nn=n+5*P*bit, aa=a+3*P*bit, bb=b+c*bit;
            Word ff=p==2 ? f+3*bit : 0, gg=p==2 ? g+2*bit : 0;
            if (aa%p + bb%p <= nn%p && !(ff & gg & 1))
                insert({nn/p,aa/p,bb/p,ff/p,gg/p});
        }
    }
    return true;
}

struct Verifier {
    int k;
    Word P;
    std::array<Word,6> counts{};

    void check(Word c, Word q) {
        if (c%15) return;
        ++counts[0];
        std::string label;
        Word kk, coefficient, value;
        require(bool(std::cin >> label >> kk >> coefficient >> value));
        require(kk == Word(k) && coefficient == c);
        if (label == "unresolved") {
            require(value == q);
            ++counts[1];
            std::cout << k << ' ' << c << ' ' << q << '\n';
        } else {
            require(label == "excluded" && (value == 2 || value == 3 || value == 5));
            require(excluded(P,c,q,value));
            ++counts[value];
        }
    }

    void lift(Word modulus, Word inverse, Word q) {
        if (modulus == 2*P) {
            check(inverse,q);
            return;
        }
        Word m = 2*modulus;
        for (Word v : {q,q+modulus}) {
            // Extend the current inverse by one digit, maintaining c*v = 1.
            Word c = inverse + ((inverse*v & modulus) ? modulus : 0);
            require((c*v) % m == 1);
            if ((3*v)%m + (2*v)%m < m) lift(m,c,v);
        }
    }

    void run() {
        lift(2,1,1);
        Word kk,total,two,three,five,remaining;
        require(bool(std::cin >> kk >> total >> two >> three >> five >> remaining));
        require(kk == Word(k) && counts[0] == total && counts[1] == remaining);
        require(counts[2] == two && counts[3] == three && counts[5] == five);
        std::cerr << "Checked " << total << " low-digit templates at " << k << std::endl;
    }
};

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    int first=std::stoi(argv[1]), last=std::stoi(argv[2]);
    if (first<3 || last>45 || first>last) return 2;
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    for (int k=first; k<=last; ++k) {
        Verifier verifier{k, Word(1)<<k};
        verifier.run();
    }
    std::string extra;
    require(!(std::cin >> extra));
}
