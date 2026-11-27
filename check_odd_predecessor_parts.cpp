// Independent census. Start with the least product of prime powers exceptional
// in a part, and split the available binary digits between the other two parts.
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using U = std::uint64_t;

U digit_sum(U n, U p) {
    U sum=0;
    while (n) {sum+=n%p; n/=p;}
    return sum;
}

bool free_in_all_bases(U n, U a, U b, U c) {
    U m=n-1;
    auto passes=[&](U p) {return digit_sum(n,p)==digit_sum(a,p)+digit_sum(b,p)+digit_sum(c,p);};
    for (U p=2; p*p<=m; ++p) {
        if (m%p) continue;
        if (!passes(p)) return false;
        while (m%p==0) m/=p;
    }
    return m==1 || passes(m);
}

int main(int argc, char** argv) {
    if (argc!=3) {
        std::cerr<<"usage: check_odd_predecessor_parts N_limit threads\n";
        return 2;
    }
    int limit=std::stoi(argv[1]), jobs=std::stoi(argv[2]);
    if (limit<3 || jobs<1) return 2;
    std::atomic<unsigned> next(2);
    std::mutex output;
    U rows=0, survivors=0;
    auto worker=[&]() {
        U count=0, keep=0;
        while (true) {
            U A=next.fetch_add(1);
            if (A*A*A>=U(limit)) break;
            for (U r=1; r<A; ++r) {
                if (std::gcd(r,A)!=1) continue;
                U x=1;
                while (r*x%A!=1) ++x;
                // h=(N-1)/A must exceed A^2 and hx is 1 modulo A.
                U h=A*A+r, a=x*h;
                for (U n=1+A*h; n<=U(limit); n+=A*A,a+=A*x) {
                    if (a & ~n) continue;
                    U remaining=n-a, m=n-1;
                    for (U b=(remaining-1)&remaining; b; b=(b-1)&remaining) {
                        U c=remaining-b;
                        if (b>=c || b%A || a*b%m || b*(b-1)%m) continue;
                        // Orthogonal idempotents a,b leave c=1-a-b modulo m.
                        if (m/std::gcd(m,b)<=A || m/std::gcd(m,c)<=A) continue;
                        ++count;
                        if (free_in_all_bases(n,a,b,c)) {
                            ++keep;
                            std::lock_guard<std::mutex> guard(output);
                            std::cout<<n<<' '<<a<<' '<<b<<' '<<c<<'\n';
                        }
                    }
                }
            }
        }
        std::lock_guard<std::mutex> guard(output);
        rows+=count;
        survivors+=keep;
    };
    std::vector<std::thread> pool;
    for (int j=0; j<jobs; ++j) pool.emplace_back(worker);
    for (auto& t:pool) t.join();
    std::cerr<<"limit "<<limit<<", odd CRT rows "<<rows<<", coprime to predecessor "<<survivors<<'\n';
}
