// Enumerate endpoint assignments by the product of one normalized part from each row.
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;
using Row = std::array<Wide,3>;

Int isqrt(Int n) {
    Int r = std::sqrt(static_cast<long double>(n));
    while ((r+1)*(r+1) <= n) ++r;
    while (r*r > n) --r;
    return r;
}

void show(Wide n) {
    if (n >= 10) show(n/10);
    std::cout << int(n%10);
}

bool avoids(Row row, Int p) {
    while (row[0] || row[1] || row[2]) {
        Int sum = 0;
        for (auto& z : row) {
            sum += Int(z%p);
            z /= p;
        }
        if (sum >= p) return false;
    }
    return true;
}

std::vector<int> primes_through(int bound) {
    std::vector<bool> prime(bound+1,true);
    std::vector<int> primes;
    for (int i=2; i<=bound; ++i) {
        if (!prime[i]) continue;
        primes.push_back(i);
        for (Int j=Int(i)*i; j<=bound; j+=i) prime[j]=false;
    }
    return primes;
}

struct Split {
    int sum;
    std::vector<int> first;
};

// Idempotents give all f modulo F with F | f(b-f) when gcd(b,F)=1.
std::vector<Split> splits(int F, int maximum, const std::vector<int>& divisors) {
    std::vector<int> residues{0};
    for (int q : divisors) {
        if (q==1) continue;
        int r=2;
        while (r*r<=q && q%r) ++r;
        if (r*r>q) r=q;
        int rest=q;
        while (rest%r==0) rest/=r;
        if (rest!=1 || (F/q)%r==0) continue;
        int other=F/q, inv=1;
        while (Int(inv)*other%q!=1) ++inv;
        int unit=other*inv, size=residues.size();
        for (int i=0; i<size; ++i) residues.push_back((residues[i]+unit)%F);
    }
    std::vector<Split> result;
    for (int b=2; b<=maximum; ++b) {
        if (b*b<4*F || std::gcd(b,F)!=1) continue;
        Split split{b,{}};
        for (int id : residues) {
            int least=Int(b)*id%F;
            for (int z=least ? least : F; z<=b/2; z+=F) split.first.push_back(z);
        }
        if (!split.first.empty()) result.push_back(std::move(split));
    }
    return result;
}

struct Search {
    int top;
    bool audit;
    std::vector<int> primes;
    std::vector<std::vector<int>> divisors;
    Int systems=0, pairs=0, good_c=0, good_pairs=0;
    Wide largest=0;
    explicit Search(int n, bool audit_pairs): top(n), audit(audit_pairs),
        // Every row scale is <= t + t*(t*t+t*t), with t<=top.
        primes(primes_through(int(isqrt(2LL*n*n*n+n))+1)), divisors(n+1) {
        for (int i=1; i<=n; ++i)
            for (int j=i; j<=n; j+=i) divisors[j].push_back(i);
    }

    bool passes(Row row, const std::array<Int,3>& assigned) const {
        for (Int number : assigned) {
            for (int r : primes) {
                if (Int(r)*r>number) break;
                if (number%r) continue;
                if (!avoids(row,r)) return false;
                do { number/=r; } while (number%r==0);
            }
            if (number>1 && !avoids(row,number)) return false;
        }
        return true;
    }

    void test(int x, int y, int F, int H, const Split& left, const Split& right,
              Int B, Int C) {
        if (std::gcd(B,C)!=1 || std::gcd(Int(F)*H,B)!=1 || std::gcd(Int(F)*H,C)!=1) return;
        Int u=y+right.sum*C, v=x+left.sum*B;
        if (std::gcd(u,v)!=1) return;
        Wide n=Wide(u)*v;
        if (n-1!=Wide(F)*H*B*C) std::abort();
        ++systems;
        largest=std::max(largest,n);
        pairs+=Int(left.first.size())*right.first.size();
        for (int j : right.first) {
            Row other{Wide(v)*y,Wide(v)*C*j,Wide(v)*C*(right.sum-j)};
            if (!passes(other,{H,C,v})) continue;
            ++good_c;
            for (int i : left.first) {
                Row first{Wide(u)*x,Wide(u)*B*i,Wide(u)*B*(left.sum-i)};
                if (audit) {
                    std::cout << "Pair";
                    for (Wide z : {Wide(x)*y, n, Wide(u), Wide(F)*B, Wide(v), Wide(H)*C}) {
                        std::cout << ' ';
                        show(z);
                    }
                    for (const auto& row : {first, other}) {
                        for (Wide z : row) {
                            std::cout << ' ';
                            show(z);
                        }
                    }
                    std::cout << '\n';
                }
                if (passes(first,{F,B,u})) {
                    ++good_pairs;
                    std::cout<<"Candidate weight "; show(n);
                    for (auto row : {first,other})
                        for (Wide z : row) {std::cout<<' ';show(z);}
                    std::cout<<std::endl;
                }
            }
        }
    }

    void roots(int x, int y, int F, int H, const Split& left, const Split& right) {
        Int b=left.sum, c=right.sum, d=F*H, e=d-b*c;
        Int target=d*x*y-e, lim=isqrt(target);
        // Take the smaller of eB-cx,eC-by. A tie belongs to the first loop.
        for (Int B=c*x/e+1; e*B-c*x<=lim; ++B) {
            Int first=e*B-c*x;
            if (target%first) continue;
            Int num=target/first+b*y;
            if (num%e==0) test(x,y,F,H,left,right,B,num/e);
        }
        if (lim*lim==target) --lim;
        for (Int C=b*y/e+1; e*C-b*y<=lim; ++C) {
            Int second=e*C-b*y;
            if (target%second) continue;
            Int num=target/second+c*x;
            if (num%e==0) test(x,y,F,H,left,right,num/e,C);
        }
    }

    void run() {
        for (int d=17; d<=top; ++d) {
            for (int F : divisors[d]) {
                int H=d/F;
                if (F>=H || std::gcd(F,H)!=1) continue;
                auto ceilroot=[](int v){return int(isqrt(4*v-1))+1;};
                auto left=splits(F,(d-1)/ceilroot(H),divisors[F]);
                if (left.empty()) continue;
                auto right=splits(H,(d-1)/left.front().sum,divisors[H]);
                for (const auto& l : left) for (const auto& r : right) {
                    if (l.sum*r.sum>=d) break;
                    for (int m=1; m<=top/d; ++m) for (int i : divisors[m])
                        roots(F*i,H*(m/i),F,H,l,r);
                }
            }
        }
    }
};

int main(int argc, char** argv) {
    int top = 10000;
    bool audit = argc == 3 && std::string_view(argv[2]) == "--audit";
    if (argc > 3 || (argc == 3 && !audit)) return 2;
    if (argc >= 2) {
        std::string_view arg(argv[1]);
        const char* end = arg.data() + arg.size();
        auto parsed = std::from_chars(arg.data(), end, top);
        if (parsed.ec != std::errc{} || parsed.ptr != end) return 2;
    }
    if (top < 1 || top > 32768) return 2;

    Search search(top, audit);
    search.run();
    std::cout << "Through " << top << ": " << search.systems << " systems, "
              << search.pairs << " row pairs, " << search.good_c << " passing second rows, "
              << search.good_pairs << " passing pairs. Largest weight ";
    show(search.largest);
    std::cout << "\n";
    return search.good_pairs ? 1 : 0;
}
