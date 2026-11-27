// Emit the unique possible row at each exponent, if its coefficient is odd.
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using Word = std::uint64_t;

unsigned halve(unsigned a, unsigned p) {
    return unsigned((Word(a) + (a % 2)*Word(p))/2);
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: check_odd_predecessor_direct P_limit Q_limit K_limit threads\n";
        return 2;
    }
    int pmax=std::stoi(argv[1]), qmax=std::stoi(argv[2]);
    int kmax=std::stoi(argv[3]), jobs=std::stoi(argv[4]);
    if (pmax<3 || qmax<pmax || kmax<1 || jobs<1) return 2;
    std::vector<bool> composite(Word(qmax)+1);
    std::vector<unsigned> primes;
    for (Word p=2; p<=Word(qmax); ++p) {
        if (composite[p]) continue;
        if (p>2) primes.push_back(unsigned(p));
        for (Word j=p*p; j<=Word(qmax); j+=p) composite[j]=true;
    }
    std::atomic<size_t> next(0);
    std::mutex output;
    Word total=0, pairs=0;
    auto worker=[&]() {
        Word local_rows=0, local_pairs=0;
        while (true) {
            size_t i=next.fetch_add(1);
            if (i>=primes.size() || primes[i]>unsigned(pmax)) break;
            unsigned p=primes[i];
            std::vector<unsigned> inverse(p);
            inverse[1]=1;
            for (unsigned t=2; t<p; ++t) inverse[t]=p-Word(p/t)*inverse[p%t]%p;
            for (size_t j=i+1; j<primes.size(); ++j) {
                unsigned q=primes[j];
                Word m=Word(p)*q;
                ++local_pairs;
                // At exponent zero, qb=1 mod p and pa=1 mod q, with pa+qb=m+1.
                unsigned b=inverse[q%p], a=q-(Word(b)*q-1)/p;
                for (Word k=1; k<=Word(kmax); ++k) {
                    a=halve(a,q);
                    b=halve(b,p);
                    Word x=Word(p)*a, y=Word(q)*b;
                    if ((y & ~m) || (x & (y | ~m))) continue;
                    ++local_rows;
                    std::lock_guard<std::mutex> lock(output);
                    std::cout<<p<<' '<<q<<' '<<k<<' '<<x<<' '<<y<<' '<<m-x-y<<'\n';
                }
            }
        }
        std::lock_guard<std::mutex> lock(output);
        total+=local_rows;
        pairs+=local_pairs;
    };
    std::vector<std::thread> pool;
    for (int j=0; j<jobs; ++j) pool.emplace_back(worker);
    for (auto& t:pool) t.join();
    std::cerr<<"pairs "<<pairs<<", weights "<<pairs*Word(kmax)<<", odd CRT rows "<<total<<'\n';
}
