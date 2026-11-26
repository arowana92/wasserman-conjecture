// Candidate endpoint assignments with a specified actual part in one row.
#include "range_factors.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

using Int = int64_t;
using Wide = __int128_t;
using Row = std::array<Int, 3>;

int t_bound(int u) {
    Int t = std::sqrt(static_cast<long double>(u)*u*u/27);
    const Wide cube = Wide(u)*u*u;
    while (27*Wide(t+1)*(t+1) <= cube) ++t;
    while (27*Wide(t)*t > cube) --t;
    return int(t);
}

int inverse(int a, int m) {
    Int x = a, y = m, u = 1, v = 0;
    while (y) {
        Int q = x/y;
        Int r = x-q*y, s = u-q*v;
        x=y; y=r; u=v; v=s;
    }
    if (x!=1) std::abort();
    return int((u%m+m)%m);
}

bool fits(Int a, Int n, int p) {
    if (p==2) return (a & n)==a;
    while (a) {
        if (a%p > n%p) return false;
        a/=p; n/=p;
    }
    return true;
}

bool avoids(const Row& row, Int n, int p) {
    return fits(row[0],n,p) && fits(row[1],n-row[0],p);
}

bool fits_all(Int a, Int n, const Factors& factors) {
    for (auto [p,w] : factors) if (!fits(a,n,p)) return false;
    return true;
}

bool avoids_all(const Row& row, Int n, const Factors& factors) {
    for (auto [p,w] : factors) if (!avoids(row,n,p)) return false;
    return true;
}

std::vector<int> divisors(const Factors& factors) {
    std::vector<int> result{1};
    for (auto [p,w] : factors) {
        size_t old = result.size();
        for (int q=w; q>1; q/=p)
            for (size_t j=0; j<old; ++j) result.push_back(q*result[j]);
    }
    return result;
}

// CRT representatives that are target at one full power and zero at the others.
std::vector<int> units(const Factors& factors, int modulus, int target) {
    std::vector<int> result;
    for (auto [p,w] : factors) {
        int other=modulus/w;
        Int k=Int(target%w)*inverse(other%w,w)%w;
        result.push_back(int(k)*other);
    }
    return result;
}

struct Search {
    FactorSieve sieve;
    int top;
    Int floor_n;
    Int systems=0, heavy_a=0, heavy_c=0, rows_a=0, rows_c=0, largest=0;

    Search(int last, Int lower) : sieve(std::max(last,t_bound(last))), top(last), floor_n(lower) {}

    void record(char side, int a, int u, Int s, int t, Int n, Row row) {
        std::sort(row.begin(),row.end());
        if (side=='A') ++rows_a; else ++rows_c;
        std::cout << side << ' ' << a << ' ' << u << ' ' << s << ' ' << t
                  << ' ' << n << ' ' << row[0] << ' ' << row[1] << ' ' << row[2] << '\n';
    }

    void companion(int a, int u, Int s, int t, Int n, const Factors& ft) {
        Int v=n/u;
        const auto basis=units(ft,t,u);
        // Distribute the atoms to three nonempty positions, in first-use order.
        auto visit=[&](auto&& self, int j, int used, Row sums) -> void {
            if (int(basis.size())-j < 3-used) return;
            if (j<int(basis.size())) {
                for (int i=0; i<std::min(3,used+1); ++i) {
                    auto next=sums;
                    next[i]+=basis[j];
                    self(self,j+1,std::max(used,i+1),next);
                }
                return;
            }
            for (Int& y : sums) y%=t;
            if (sums[0]+sums[1]+sums[2]!=u) return;
            for (Int& y : sums) y*=v;
            if (avoids_all(sums,n,ft)) record('C',a,u,s,t,n,sums);
        };
        visit(visit,0,0,{});
    }

    void anchored(int a, int u, int f, int b, int t, Int n,
                  const Factors& fu, const Factors& ff, const Factors& fb) {
        const Int step=Int(u)*b, sum=(n-a)/step;
        // Since S>V, sum<F. Every possible first summand is a least residue.
        if (sum>=f) std::abort();
        auto basis=units(ff,f,int(sum));
        std::vector<int> residues{0};
        for (int delta : basis) {
            size_t old=residues.size();
            for (size_t i=0; i<old; ++i) residues.push_back((residues[i]+delta)%f);
        }
        for (int h : residues) {
            if (!h || 2*h>sum) continue;
            Row row{a,step*h,step*(sum-h)};
            if (avoids_all(row,n,fu) && avoids_all(row,n,ff) && avoids_all(row,n,fb))
                record('A',a,u,Int(f)*b,t,n,row);
        }
    }

    void run(int first) {
        for (int a=first; a<=top; ++a) {
            auto da=divisors(sieve.factor(a));
            auto pred=divisors(sieve.factor(a-1));
            for (int u : da) {
                if (u<3) continue;
                int x=a/u, hi=t_bound(u);
                Factors fu=sieve.factor(u);
                for (int f : divisors(sieve.factor(x))) {
                    if (std::gcd(f,u)!=1) continue;
                    int inv=inverse(f%u,u);
                    Factors ff=sieve.factor(f);
                    for (int b : pred) {
                        Int s=Int(f)*b;
                        if (s*hi<floor_n || std::gcd(f,b)!=1) continue;
                        int k=(a-1)/b;
                        Int lo=std::max({Int(1),floor_n/s,(Int(k)+2*u+f-1)/f});
                        if (lo>hi) continue;
                        Factors fb=sieve.factor(b);
                        Int rem=Int(k)*inv%u;
                        for (Int t=lo+(rem-lo%u+u)%u; t<=hi; t+=u) {
                            Int n=s*t+1, v=n/u;
                            if (n<=floor_n || std::gcd(s,t)!=1 || std::gcd(Int(u),v)!=1) continue;
                            ++systems; largest=std::max(largest,n);
                            if (!fits_all(a,n,fu) || !fits_all(a,n,ff) || !fits_all(a,n,fb)) continue;
                            Factors ft=sieve.factor(int(t));
                            bool feasible=true;
                            for (auto [p,w] : ft) if (2LL*w+1>u) feasible=false;
                            if (!feasible) continue;
                            if (t>u) {
                                if (ft.size<3) continue;
                                ++heavy_c;
                                companion(a,u,s,int(t),n,ft);
                            } else {
                                ++heavy_a;
                                anchored(a,u,f,b,int(t),n,fu,ff,fb);
                            }
                        }
                    }
                }
            }
        }
        std::cerr << "range " << first << ' ' << top << " above " << floor_n
                  << " systems " << systems << " largest " << largest << " heavy A/C "
                  << heavy_a << ' ' << heavy_c << " passing rows A/C " << rows_a << ' ' << rows_c << '\n';
    }
};

int main(int argc, char** argv) {
    if (argc==2 && std::string(argv[1])=="--example") {
        Search search(283,0);
        FactorSieve table(887603);
        const int n=222788354, s=887603, u=283;
        for (int a : {86617810,45491401,90679143}) {
            int f=std::gcd(a/u,s), b=s/f;
            search.anchored(a,u,f,b,251,n,table.factor(u),table.factor(f),table.factor(b));
        }
        search.companion(787238,787238,251,s,n,table.factor(s));
        return 0;
    }
    if (argc!=4) return 2;
    int first=std::stoi(argv[1]), last=std::stoi(argv[2]);
    Int lower=std::stoll(argv[3]);
    if (first<3 || last<first || last>1000000 || lower<0) return 2;
    Search search(last,lower);
    search.run(first);
}
