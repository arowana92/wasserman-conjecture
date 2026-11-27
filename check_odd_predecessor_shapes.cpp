// Generate every binary allocation with a possible odd weight prime at N=1+2^k*p*q.
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Word = std::uint64_t;

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: check_odd_predecessor_shapes P_limit Q_limit threads\n";
        return 2;
    }
    int p_limit = std::stoi(argv[1]), q_limit = std::stoi(argv[2]), threads = std::stoi(argv[3]);
    if (p_limit < 3 || q_limit < p_limit || threads < 1) return 2;
    std::vector<bool> prime(Word(q_limit)+1, true);
    std::vector<int> odd_primes;
    for (Word p=2; p<=Word(q_limit); ++p) {
        if (!prime[p]) continue;
        if (p>2) odd_primes.push_back(p);
        for (Word j=Word(p)*p; j<=Word(q_limit); j+=p) prime[j]=false;
    }
    std::atomic<size_t> next(0);
    std::mutex output;
    Word pairs=0, candidates=0, candidate_pairs=0;
    auto worker = [&]() {
        Word local_pairs=0, local_rows=0, local_positive=0;
        while (true) {
            size_t ip=next.fetch_add(1);
            if (ip>=odd_primes.size() || odd_primes[ip]>p_limit) break;
            int p=odd_primes[ip];
            std::vector<bool> powers(p);
            for (int u=1; !powers[u]; u=(2*Word(u))%p) powers[u]=true;
            for (size_t iq=ip+1; iq<odd_primes.size(); ++iq) {
                int q=odd_primes[iq];
                Word m=Word(p)*q;
                ++local_pairs;
                bool found=false;
                for (int b=3; b<p; ++b) {
                    int odd_b=b;
                    while (odd_b%2==0) odd_b/=2;
                    Word y=Word(b)*q;
                    if (odd_b==1 || (y & ~m) || !powers[y%p]) continue;
                    Word available=m-y;
                    for (Word x=(available-1)&available; x; x=(x-1)&available) {
                        if (x%p) continue;
                        Word a=x/p;
                        if (std::gcd(a,Word(odd_b))==1) continue;
                        found=true;
                        ++local_rows;
                        std::lock_guard<std::mutex> lock(output);
                        std::cout<<p<<' '<<q<<' '<<a<<' '<<b<<'\n';
                    }
                }
                local_positive+=found;
            }
        }
        std::lock_guard<std::mutex> lock(output);
        pairs+=local_pairs;
        candidates+=local_rows;
        candidate_pairs+=local_positive;
    };
    std::vector<std::thread> pool;
    for (int j=0; j<threads; ++j) pool.emplace_back(worker);
    for (auto& t:pool) t.join();
    std::cerr<<"pairs "<<pairs<<", pairs passing binary tests "<<candidate_pairs
             <<", rows "<<candidates<<'\n';
}
