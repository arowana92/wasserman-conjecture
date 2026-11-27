// Exhaust both companion families for 10P^2+7P+1 by CRT at 10P+7.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using Int = std::int64_t;
struct Counts {
    Int parameters = 0, allocations = 0, rows = 0, parity_pass = 0;
};

Int inverse(Int a, Int modulus) {
    Int b = modulus, x = 1, y = 0;
    while (b != 0) {
        Int q = a/b, r = a%b, next = x-q*y;
        a = b; b = r; x = y; y = next;
    }
    if (a != 1) throw std::logic_error("inverse requires a unit");
    return (x%modulus+modulus)%modulus;
}

bool prime(int n) {
    for (int p = 2; p*p <= n; ++p) if (n%p == 0) return false;
    return n >= 2;
}

void check_orientation(int P, const std::vector<int>& factors, int t, Counts& counts) {
    const Int T = 1+Int(t)*P, W = 1+Int(7-t)*P, H = 10LL*P+7;
    if (factors.size() < 3) return;
    ++counts.parameters;
    std::vector<Int> terms;
    for (int power : factors) {
        Int rest = H/power;
        terms.push_back(rest*(W*inverse(rest, power)%power));
    }
    auto visit = [&](auto&& self, size_t index, int used, std::array<Int,3> f) -> void {
        if (used+terms.size()-index < 3) return;
        if (index < terms.size()) {
            for (int i = 0; i < std::min(3,used+1); ++i) {
                auto next = f;
                next[i] = (next[i]+terms[index])%H;
                self(self,index+1,std::max(used,i+1),next);
            }
            return;
        }
        ++counts.allocations;
        if (f[0]+f[1]+f[2] != W) return;
        ++counts.rows;
        const std::array<Int,3> C = {T*f[0], T*f[1], T*f[2]};
        bool need_odd = t == 2 || P%2 == 1;
        if (need_odd && ((C[0]&C[1]) || (C[0]&C[2]) || (C[1]&C[2]))) return;
        ++counts.parity_pass;
        std::cout << "Unexcluded P=" << P << ", t=" << t << ", companion "
                  << C[0] << " " << C[1] << " " << C[2] << "\n";
    };
    visit(visit,0,0,{0,0,0});
}

int main(int argc, char** argv) {
    const int limit = argc == 2 ? std::stoi(argv[1]) : 10000000;
    if (limit < 1 || limit > 100000000) throw std::invalid_argument("bound must be in [1,100000000]");
    std::vector<std::pair<int,int>> sieve;
    for (int p = 3; Int(p)*p <= 10LL*limit+7; p += 2) {
        if (p == 5 || !prime(p)) continue;
        sieve.push_back({p, int((p-7%p)*inverse(10,p)%p)});
    }
    std::array<Counts,2> counts;
    constexpr int block = 65536;
    for (int start = 1; start <= limit; start += block) {
        int length = std::min(block,limit-start+1);
        std::vector<int> left(length);
        std::vector<std::vector<int>> factors(length);
        for (int i = 0; i < length; ++i) left[i] = 10*(start+i)+7;
        for (auto [p,residue] : sieve) {
            for (int i = (residue-start%p+p)%p; i < length; i += p) {
                int power = 1;
                do { power *= p; left[i] /= p; } while (left[i]%p == 0);
                factors[i].push_back(power);
            }
        }
        for (int i = 0; i < length; ++i) {
            if (left[i] > 1) factors[i].push_back(left[i]);
            int P = start+i;
            check_orientation(P,factors[i],2,counts[0]);
            if (P%2 || (P & (P-1)) == 0)
                check_orientation(P,factors[i],5,counts[1]);
        }
    }
    for (size_t i = 0; i < counts.size(); ++i) {
        const auto& c = counts[i];
        std::cout << "P through " << limit << ", t=" << (i == 0 ? 2 : 5)
                  << ", parameters with three available bases: " << c.parameters
                  << ", partitions: " << c.allocations << ", rows: " << c.rows
                  << ", after parity: " << c.parity_pass << "\n";
    }
    return counts[0].parity_pass != 0 || counts[1].parity_pass != 0;
}
