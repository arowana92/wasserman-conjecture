// Search endpoint assignments with a prescribed bound on the spread of a nonconstant row.
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

int upper_t(int u) {
    const Int cube = Int(u)*u*u;
    int t = int(std::sqrt(static_cast<long double>(cube)/27));
    while (27*Int(t+1)*(t+1) <= cube) ++t;
    while (27*Int(t)*t > cube) --t;
    return t;
}

int inverse_mod(Int n, int modulus) {
    Int a=n%modulus, b=modulus, x=1, y=0;
    while (b) {
        Int q=a/b, r=a%b, z=x-q*y;
        a=b; b=r; x=y; y=z;
    }
    if (a != 1) std::abort();
    return int((x%modulus+modulus)%modulus);
}

bool avoids(std::array<Int,3> row, int p) {
    while (row[0] || row[1] || row[2]) {
        int digits=0;
        for (Int& part : row) {
            digits += int(part%p);
            part /= p;
        }
        if (digits >= p) return false;
    }
    return true;
}

int odd_part(int n) {
    return n/(n & -n);
}

struct Search {
    int bound;
    Int floor;
    FactorSieve sieve;
    std::vector<std::vector<int>> odd_divisors;
    Int systems[2]{}, passing[2]{}, largest[2]{};
    bool emit;

    Search(int limit, Int lower, bool output) : bound(limit), floor(lower), sieve(limit+1),
            odd_divisors(limit+2), emit(output) {
        for (int i=1; i<=limit+1; i+=2)
            for (int j=i; j<=limit+1; j+=i) odd_divisors[j].push_back(i);
    }

    bool first_row_passes(const std::array<Int,3>& row, const std::array<int,3>& bins,
                          const Factors& weight_factors, const std::array<Factors,3>& bin_factors) {
        for (const Atom& atom : weight_factors)
            if (!avoids(row,atom.prime)) return false;
        for (int i=0; i<3; ++i) {
            if ((row[i]-1)%bins[i] || row[(i+1)%3]%bins[i] || row[(i+2)%3]%bins[i])
                return false;
            for (const Atom& atom : bin_factors[i])
                if (!avoids(row,atom.prime)) return false;
        }
        return true;
    }

    void progression(int d, int f, int u, std::array<int,3> bins, Int s, int hi) {
        int rem = (u-inverse_mod(s,u))%u;
        Int lo = (floor-1)/s+1;
        const Factors weight_factors=sieve.factor(u);
        const std::array<Factors,3> bin_factors{sieve.factor(bins[0]),sieve.factor(bins[1]),sieve.factor(bins[2])};
        for (Int t=lo+(rem-lo%u+u)%u; t<=hi; t+=u) {
            Int n=s*t+1, v=n/u, diff=n-2*d-f;
            if (diff%(3*u) || std::gcd(v,Int(u))!=1 || std::gcd(s,t)!=1) continue;
            Int a=diff/3;
            std::array<Int,3> row{a,a+d,a+d+f};
            bool distinct=d && f;
            if (a<=0) std::abort(); // n > 4*bound
            ++systems[distinct];
            largest[distinct]=std::max(largest[distinct],n);
            if (!first_row_passes(row,bins,weight_factors,bin_factors)) continue;
            ++passing[distinct];
            if (emit) std::cout<<n<<' '<<u<<' '<<s<<' '<<t<<' '<<a<<' '<<d<<' '<<f<<'\n';
        }
    }

    void gaps(int g, int d, int f) {
        int b1=odd_part(std::gcd(f,d+1)), b2=odd_part(std::gcd(d+f,d-1));
        int b3=odd_part(std::gcd(d,f-1));
        Int capacity=Int(b1)*b2*b3;
        if (capacity*upper_t(g) < floor) return;
        for (int u : odd_divisors[g]) {
            if (u<7 || std::gcd(capacity,Int(u))!=1) continue;
            int hi=upper_t(u);
            if (capacity*hi < floor) continue;
            for (int i : odd_divisors[b1])
                for (int j : odd_divisors[b2])
                    for (int k : odd_divisors[b3]) {
                        Int s=Int(i)*j*k;
                        if (s*hi >= floor) progression(d,f,u,{i,j,k},s,hi);
                    }
        }
    }

    // Choose U and the three assigned bins from the span, independently of the
    // primitive gap pair. The normalized first gap then lies in one CRT class.
    void by_span() {
        for (int w=7; w<=bound; ++w) {
            for (int u : odd_divisors[w]) {
                if (u<7) continue;
                int hi=upper_t(u), h=w/u;
                for (int left : odd_divisors[w+1]) {
                    for (int right : odd_divisors[w-1]) {
                        Int outer=Int(left)*right;
                        if (outer>w+1 || outer*hi*h<floor) continue;
                        int initial=(left-inverse_mod(Int(u)*right,left))%left;
                        Int x0=Int(right)*initial;
                        if (x0>h) continue;
                        for (int middle : odd_divisors[h]) {
                            Int s=outer*middle;
                            if (s*hi<floor || std::gcd(u,middle)!=1) continue;
                            // u*(x0+outer*k) = 1 modulo the middle bin.
                            Int k=(1-Int(u)*x0)%middle;
                            k=k*inverse_mod(Int(u)*outer,middle)%middle;
                            if (k<0) k+=middle;
                            for (Int x=x0+outer*k; x<=h; x+=s) {
                                int d=int(u*x);
                                progression(d,w-d,u,{left,middle,right},s,hi);
                            }
                        }
                    }
                }
            }
        }
    }

    void by_primitive_gaps() {
        for (int g=7; g<=bound; ++g) {
            if (odd_part(g)<7) continue;
            gaps(g,0,g);
            gaps(g,g,0);
            int hi=upper_t(odd_part(g));
            // B1*B2, B2*B3 <= 2W+1 and B1*B3 <= W+1.
            if (Wide(2*bound+1)*(2*bound+1)*(bound+1)*hi*hi < Wide(floor)*floor) continue;
            int threshold=int(std::cbrt(static_cast<long double>(floor)/hi));
            while (Int(threshold)*threshold*threshold*hi<floor) ++threshold;
            while (threshold>1 && Int(threshold-1)*(threshold-1)*(threshold-1)*hi>=floor) --threshold;
            for (int h=2; h<=bound/g; ++h) {
                if (Int(h)*h*h*hi < 4*floor) continue;
                int w=g*h;
                // One of the three maximal odd bins must reach the cube-root threshold.
                std::array<int,3> moduli{w+1,h,w-1};
                for (int i=0; i<3; ++i) {
                    const auto& ds=odd_divisors[moduli[i]];
                    for (auto it=std::lower_bound(ds.begin(),ds.end(),threshold);
                         it!=ds.end() && *it<=h; ++it) {
                        int bin=*it, start=0;
                        if (i==0) start=h%bin;
                        if (i==1) {
                            if (std::gcd(g,bin)!=1) continue;
                            start=inverse_mod(g,bin);
                        }
                        for (int r=start?start:bin; r<h; r+=bin) {
                            if (std::gcd(r,h)!=1) continue;
                            int d=g*r, f=w-d;
                            std::array<int,3> bs{odd_part(std::gcd(f,d+1)),
                                odd_part(std::gcd(w,d-1)),odd_part(std::gcd(d,f-1))};
                            auto largest_bin=std::max_element(bs.begin(),bs.end());
                            if (largest_bin==bs.begin()+i && *largest_bin==bin) gaps(g,d,f);
                        }
                    }
                }
            }
        }
    }
};

int main(int argc,char** argv) {
    int bound=argc>1?std::stoi(argv[1]):1000000;
    Int floor=argc>2?std::stoll(argv[2]):1000000000LL;
    if (bound<7 || bound>1000000 || floor<4LL*bound || floor>1000000000000LL) return 2;
    bool emit=false, by_span=false;
    for (int i=3;i<argc;++i) {
        if (std::string(argv[i])=="--emit") emit=true;
        else if (std::string(argv[i])=="--by-span") by_span=true;
        else return 2;
    }
    Search search(bound,floor,emit);
    if (by_span) search.by_span(); else search.by_primitive_gaps();
    for (int i=0;i<2;++i)
        std::cerr<<(i?"distinct":"repeated")<<" systems "<<search.systems[i]<<" passing "
                 <<search.passing[i]<<" largest "<<search.largest[i]<<'\n';
}
