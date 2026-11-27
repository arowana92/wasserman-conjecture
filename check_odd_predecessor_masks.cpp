// Reconstruct the two unshifted slopes from the binary mask of their sum.
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using U64 = std::uint64_t;

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    int pmax=std::stoi(argv[1]), qmax=std::stoi(argv[2]), threads=std::stoi(argv[3]);
    if (pmax<3 || qmax<pmax || threads<1) return 2;
    std::vector<int> primes;
    for (U64 n=3; n<=U64(qmax); n+=2) {
        bool prime=true;
        for (auto f:primes) {
            if (f>n/f) break;
            if (n%f==0) {prime=false; break;}
        }
        if (prime) primes.push_back(n);
    }
    std::atomic<size_t> next_q(0);
    std::mutex output;
    U64 pairs=0, rows=0;
    auto worker=[&]() {
        U64 local_pairs=0, local_rows=0;
        while (true) {
            size_t at=next_q.fetch_add(1);
            if (at>=primes.size()) break;
            int q=primes[at];
            for (int p:primes) {
                if (p>pmax || p>=q) break;
                ++local_pairs;
                U64 m=U64(p)*q;
                // t=x+y is 2^{-k} at both odd primes. Split its CRT residues
                // by y=0 mod q and y=t mod p. The inverse only depends on q mod p.
                int inv=1;
                while (U64(inv)*q%p!=1) ++inv;
                std::vector<bool> orbit(p);
                int v=1;
                do {orbit[v]=true; v=(v&1)?(U64(v)+p)/2:v/2;} while (v!=1);
                for (U64 t=(m-1)&m; t; t=(t-1)&m) {
                    U64 b=(t%p)*inv%p, y=b*q;
                    if (!b || y>=t || (y & ~t) || !orbit[t%p]) continue;
                    U64 a=(t-y)/p;
                    U64 factor=std::gcd(a,b);
                    if ((factor & (factor-1))==0) continue;
                    ++local_rows;
                    std::lock_guard<std::mutex> guard(output);
                    std::cout<<p<<' '<<q<<' '<<a<<' '<<b<<'\n';
                }
            }
        }
        std::lock_guard<std::mutex> guard(output);
        pairs+=local_pairs; rows+=local_rows;
    };
    std::vector<std::thread> pool;
    for (int i=0; i<threads; ++i) pool.emplace_back(worker);
    for (auto& t:pool) t.join();
    std::cerr<<"pairs "<<pairs<<", binary shapes "<<rows<<'\n';
}
