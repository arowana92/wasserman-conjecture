// Direct orbit enumeration and an independent carry test using two parts and the total.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using I = std::int64_t;
using U = std::uint64_t;

I power_mod(I a, I n, I modulus) {
    I answer=1;
    while (n) {
        if (n%2) answer = (__int128(answer)*a)%modulus;
        a = (__int128(a)*a)%modulus;
        n/=2;
    }
    return answer;
}

bool fits(I x, I y, I n, int base) {
    while (x || y || n) {
        if (x%base + y%base > n%base) return false;
        x/=base; y/=base; n/=base;
    }
    return true;
}

bool possible(int q, int z, int p, int k, int stride) {
    // At 3 and 5, v_p(2^stride-1)=1+v_p(stride) for 4|stride.
    I modulus=p;
    for (int j=stride; j%p==0; j/=p) modulus*=p;
    if (p==3 && modulus<9) modulus=9;

    const I first=3*q, total=5*q;
    // For a low input t with denominator D, x=floor(3qt/D).
    // Given x, y differs from floor(zx/(3q)) by 0 or 1 and
    // n differs from floor(5x/3) by 0, 1 or 2 since D>=3.
    std::vector<U> seen(((first+1)*6+63)/64);
    std::vector<U> queue;
    auto add = [&](I x,I y,I n) {
        I dy=y-z*x/first, dn=n-5*x/3;
        if (dy<0 || dy>1 || dn<0 || dn>2) throw std::logic_error("carry state");
        U id=6*x+3*dy+dn;
        U mask=U(1)<<(id%64);
        if (!(seen.at(id/64)&mask)) {
            seen[id/64]|=mask;
            queue.push_back(id);
        }
    };
    I start=power_mod(2,k,modulus), ratio=power_mod(2,stride,modulus);
    I w=start;
    do {
        I x=first*w, y=I(z)*w, n=total*w+1;
        if ((p!=3 || n%9) && fits(x%modulus,y%modulus,n%modulus,p)) {
            add(x/modulus,y/modulus,n/modulus);
        }
        w=w*ratio%modulus;
    } while (w!=start);
    for (std::size_t i=0; i<queue.size(); ++i) {
        I index=queue[i], x=index/6;
        I y=z*x/first+(index%6)/3, n=5*x/3+index%3;
        if (fits(x,y,n,p)) return true;
        for (int digit=0; digit<p; ++digit) {
            I xx=x+first*digit, yy=y+I(z)*digit, nn=n+total*digit;
            if (xx%p+yy%p > nn%p) continue;
            add(xx/p,yy/p,nn/p);
        }
    }
    return false;
}

void search(int bits) {
    int bound=1<<bits;
    std::vector<int> factor(bound);
    for (int p=2; p<bound; ++p) {
        if (factor[p]) continue;
        for (int j=p; j<bound; j+=p) if (!factor[j]) factor[j]=p;
    }
    I powers=0, classes=0, reject3=0, reject5=0;
    for (int q=11; q<bound; q+=8) {
        if ((3*q & 2*q) || __builtin_popcount(unsigned(q))<7) continue;
        int p=factor[q], t=q;
        while (t%p==0) t/=p;
        if (t!=1 || p<=5) continue;
        int first=0;
        while (first<4 && (I(q)<<first)%15!=7) ++first;
        if (first==4) continue;
        ++powers;
        int half=(q+1)/2;
        I start=power_mod(half,first,q), inverse16=power_mod(half,4,q);
        I a=start;
        int k=first;
        std::vector<std::pair<int,int>> rows;
        do {
            int z=a%2 ? a+q : a;
            if (z%15==0 && !(z & ~(2*q))) {
                // Recover 2^k directly and check every floor at powers
                // p^(e+1),...,p^(2e) for the distinguished 3q cut.
                I w=power_mod(2,k,q);
                bool pass=true;
                for (I g=p; g<=q; g*=p)
                    pass &= 3*w/g+2*w/g == 5*w/g;
                if (pass) rows.emplace_back(z,k);
            }
            a=a*inverse16%q; k+=4;
        } while (a!=start);
        int stride=k-first;
        for (auto [z,r]:rows) {
            ++classes;
            if (!possible(q,z,3,r,stride)) ++reject3;
            else if (!possible(q,z,5,r,stride)) ++reject5;
            else std::cout << "survivor " << q << ' ' << z << ' ' << r << ' ' << stride << '\n';
        }
    }
    std::cout << "prime powers, classes, excluded at 3, excluded at 5\n"
              << powers << ' ' << classes << ' ' << reject3 << ' ' << reject5 << '\n';
}

int main(int argc,char** argv) {
    int bits=argc==2 ? std::stoi(argv[1]):21;
    if (argc>2 || bits<3 || bits>24) return 2;
    search(bits);
}
