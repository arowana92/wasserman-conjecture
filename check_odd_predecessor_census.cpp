// Enumerate the orthogonal idempotents of N-1 that fit the binary digits of N.
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using Int = std::int64_t;

int inverse(int a, int m) {
    Int x=a, y=m, u=1, v=0;
    while (y) {
        Int k=x/y;
        x-=k*y; std::swap(x,y);
        u-=k*v; std::swap(u,v);
    }
    return int((u%m+m)%m);
}

bool factorial_floors(int n, int a, int b, int c, int p) {
    for (Int w=p; w<=n; w*=p)
        if (n/w != a/w+b/w+c/w) return false;
    return true;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr<<"usage: check_odd_predecessor_census N_limit threads\n";
        return 2;
    }
    int limit=std::stoi(argv[1]), jobs=std::stoi(argv[2]);
    if (limit<3 || jobs<1) return 2;
    // An odd composite below INT_MAX has a smallest factor below 46341.
    std::vector<std::uint16_t> smallest(limit/2+1);
    for (int p=3; Int(p)*p<=limit; p+=2)
        if (!smallest[p/2])
            for (Int j=Int(p)*p; j<=limit; j+=2*p)
                if (!smallest[j/2]) smallest[j/2]=std::uint16_t(p);

    std::atomic<Int> next(3);
    std::mutex output;
    Int odd_rows=0, surviving=0;
    auto worker=[&]() {
        Int count=0, keep=0;
        // The product of the first ten primes exceeds INT_MAX.
        std::array<int,9> bases{}, powers{};
        std::array<int,512> residues{}, masks{};
        while (true) {
            Int start=next.fetch_add(32768);
            if (start>limit) break;
            for (Int nn=start; nn<std::min(start+32768,Int(limit)+1); ++nn) {
                int n=int(nn), m=n-1, atoms=0, rest=m;
                while (rest>1) {
                    int p=(rest%2==0) ? 2 : smallest[rest/2];
                    if (!p) p=rest;
                    int w=1;
                    do {rest/=p; w*=p;} while (rest%p==0);
                    bases[atoms]=p;
                    powers[atoms++]=w;
                }
                if (atoms<3) continue;
                int size=1;
                for (int i=0; i<atoms; ++i) {
                    int w=powers[i], scale=m/w;
                    int e=scale*inverse(scale%w,w);
                    for (int j=0; j<size; ++j)
                        residues[j+size]=(residues[j]>=m-e) ? residues[j]-(m-e) : residues[j]+e;
                    size*=2;
                }
                int possible=0;
                for (int mask=1; mask<size-1; ++mask)
                    if (!(residues[mask] & ~n)) masks[possible++]=mask;
                // The three part masks and their three pairwise unions all fit.
                if (possible<6) continue;
                for (int i=0; i<possible; ++i) for (int j=i+1; j<possible; ++j) {
                    int A=masks[i], B=masks[j], C=(size-1)^(A|B);
                    if ((A&B) || C<=B) continue;
                    int a=residues[A], b=residues[B], c=residues[C];
                    if ((a&b) || Int(a)+b+c != n) continue;
                    ++count;
                    bool free=true;
                    for (int t=0; t<atoms && free; ++t)
                        free=factorial_floors(n,a,b,c,bases[t]);
                    if (free) {
                        ++keep;
                        std::lock_guard<std::mutex> guard(output);
                        std::cout<<n<<' '<<a<<' '<<b<<' '<<c<<'\n';
                    }
                }
            }
        }
        std::lock_guard<std::mutex> guard(output);
        odd_rows+=count;
        surviving+=keep;
    };
    std::vector<std::thread> pool;
    for (int j=0; j<jobs; ++j) pool.emplace_back(worker);
    for (auto& t:pool) t.join();
    std::cerr<<"limit "<<limit<<", odd CRT rows "<<odd_rows<<", coprime to predecessor "<<surviving<<'\n';
}
