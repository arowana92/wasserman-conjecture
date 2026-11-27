// Enumerate all lowest-layer pairs having an odd coefficient at their candidate weight.
#include "range_factors.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

using Int = int64_t;
using Row = std::array<int, 3>;
using Record = std::tuple<Row, Row, Int>;

int inverse(int a, int b) {
    Int r = a, s = b, x = 1, y = 0;
    while (s) {
        Int q = r/s;
        std::tie(r, s) = std::make_pair(s, r-q*s);
        std::tie(x, y) = std::make_pair(y, x-q*y);
    }
    if (r != 1) throw std::logic_error("nonunit in CRT");
    return int((x % b + b) % b);
}

int group_bound(int limit) {
    const Int cube = Int(limit)*limit*limit;
    int t = int(std::sqrt(static_cast<double>(cube/27)));
    while (Int(t+1)*(t+1) <= cube/27) ++t;
    while (Int(t)*t > cube/27) --t;
    return t;
}

bool odd_row(const Row& x, Int scale) {
    const Int a = scale*x[0], b = scale*x[1], c = scale*x[2];
    return !(a & b) && !((a | b) & c);
}

struct Search {
    int limit, max_t;
    FactorSieve sieve;
    std::vector<int> trial_primes;
    std::vector<std::vector<Row>> small_residues;
    std::set<Record> accepted;

    explicit Search(int x) : limit(x), max_t(group_bound(x)),
        sieve(std::max(x,max_t+1)), small_residues(x+1) {
        for (int p=2; Int(p)*p<=Int(x)*x/27+1; ++p)
            if (sieve.prime(p)) trial_primes.push_back(p);
    }

    Factors scale_factors(int m) const {
        if (m<=max_t+1) return sieve.factor(m);
        Factors result;
        for (int p : trial_primes) {
            if (Int(p)*p>m) break;
            if (m%p) continue;
            int w=1;
            do { m/=p; w*=p; } while (m%p==0);
            result.atoms[result.size++] = {p,w};
        }
        if (m>1) result.atoms[result.size++] = {m,m};
        return result;
    }

    std::vector<Row> make_residues(int t) const {
        // Form unlabeled partitions of its full prime powers into at most three groups.
        // Entries are CRT idempotents with target one on the indicated group.
        std::vector<Row> groups{{0,0,0}};
        for (auto [p, w] : sieve.factor(t)) {
            int unit = (t/w)*inverse((t/w) % w, w);
            std::vector<Row> next;
            for (Row old : groups) {
                for (int i=0; i<3; ++i) {
                    Row parts = old;
                    parts[i] = int((Int(old[i])+unit) % t);
                    next.push_back(parts);
                    if (old[i] == 0) break;
                }
            }
            groups.swap(next);
        }
        return groups;
    }

    void companions(const Row& x, Int n, bool odd_x, int v, int t,
                    const std::vector<Row>& classes) {
        const Int scale = n/v;
        for (Row r : classes) {
            for (int& a : r) {
                a = int(Int(a)*v % t);
                if (a == 0) a = t;
            }
            for (int a=r[0]; a+r[1]+r[2]<=v; a+=t) {
                const Int part = a*scale;
                if (!odd_x && (part & n) != part) continue;
                for (int b=r[1]; a+b+r[2]<=v; b+=t) {
                    Row y{a, b, v-a-b};
                    if (!odd_x && !odd_row(y, scale)) continue;
                    if (std::gcd(std::gcd(a, b), v) != 1) continue;
                    Row first = x;
                    std::sort(first.begin(), first.end());
                    std::sort(y.begin(), y.end());
                    if (y < first) std::swap(first, y);
                    accepted.emplace(first, y, n);
                }
            }
        }
    }

    void process(const Row& x, int u, int s) {
        // If T >= v, both rows qualify as the first row. Enumerate that case
        // only when the first row has the odd coefficient.
        const int low_a = x[0] & -x[0], low_b = x[1] & -x[1], low_c = x[2] & -x[2];
        const bool distinct = low_a != low_b && low_a != low_c && low_b != low_c;
        const int upper = distinct ? max_t : std::min(max_t, limit-1);
        for (int t=u-inverse(s%u, u); t<=upper; t+=u) {
            const Int n = Int(s)*t+1;
            const int m = int(n/u);
            const bool odd_x = distinct && odd_row(x, m);
            if ((!odd_x && t>=limit) || std::gcd(s, t)!=1) continue;
            std::vector<int> choices{1};
            for (auto [p,w] : scale_factors(m)) {
                size_t count = choices.size();
                for (int power=w; power>1; power/=p)
                    for (size_t i=0; i<count; ++i)
                        if (choices[i] <= limit/power) choices.push_back(choices[i]*power);
            }
            std::vector<Row> large_classes;
            for (int v : choices) {
                if (v<3 || (!odd_x && (!(v%2) || t>=v)) ||
                    Int(t)*t > Int(v)*v*v/27 || std::gcd(v,u)!=1) continue;
                auto& classes = t<=limit ? small_residues[t] : large_classes;
                if (classes.empty()) classes = make_residues(t);
                companions(x, n, odd_x, v, t, classes);
            }
        }
    }

    void run(int first_d, int last_d) {
        for (int d=first_d; d<=last_d && Int(3)*d*d<=limit; ++d) {
            for (int e=d+1; e*e+2*d*e<=limit; ++e) {
                if (std::gcd(d,e)!=1) continue;
                for (int f=e+1; d*e+(d+e)*f<=limit; ++f) {
                    if (std::gcd(d*e,f)!=1) continue;
                    const int a=e*f, b=d*f, c=d*e, s=d*e*f;
                    const int cap=std::min(s,limit);
                    for (int i=a; i+b+c<=cap; i+=a)
                        for (int j=b; i+j+c<=cap; j+=b) {
                            int h=std::gcd(i,j);
                            for (int k=c; i+j+k<=cap; k+=c)
                                if (std::gcd(h,k)==1) process({i,j,k},i+j+k,s);
                        }
                }
            }
        }
        for (const auto& [x,y,n] : accepted)
            std::cout << n << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << ' '
                      << y[0] << ' ' << y[1] << ' ' << y[2] << '\n';
    }
};

int main(int argc, char** argv) {
    if (argc!=2 && argc!=4) return 2;
    int limit=std::stoi(argv[1]);
    // Keep m <= X^2/27 + 1 below 10^9 for the integer factor representation.
    if (limit<3 || limit>160000) return 2;
    int first=argc==4 ? std::stoi(argv[2]) : 2;
    int last=argc==4 ? std::stoi(argv[3]) : limit;
    if (first<2 || first>last || first>limit) return 2;
    Search search(limit);
    search.run(first,last);
}
