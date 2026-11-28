// Enumerate primitive triples by coordinates and compute the available prime powers.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <unordered_set>
#include <vector>

using Row = std::array<int,3>;
using Powers = std::vector<std::pair<int,int>>;

Powers powers(int h) {
    Powers result;
    for (int p=2; p*p<=h; ++p) {
        int q=1;
        while (h%p==0) { h/=p; q*=p; }
        if (q>1) result.emplace_back(p,q);
    }
    if (h>1) result.emplace_back(h,h);
    return result;
}

bool carry_free(Row row,int p) {
    while (row[0] || row[1] || row[2]) {
        int sum=0;
        for (int &x:row) { sum+=x%p; x/=p; }
        if (sum>=p) return false;
    }
    return true;
}

bool possible(const Row &row,int p,int q) {
    const int u=row[0]+row[1]+row[2];
    int r=1;
    while ((1LL*r*u)%q!=1) ++r;
    Row first{};
    for (int i=0;i<3;++i) first[i]=row[i]*r;
    Row suffix=first;
    for (int &x:suffix) x%=q;
    if (!carry_free(suffix,p)) return false;
    for (int &x:first) x/=q;
    std::vector<Row> todo{first};
    std::unordered_set<long long> seen;
    for (size_t k=0;k<todo.size();++k) {
        const Row state=todo[k];
        if (carry_free(state,p)) return true;
        const int radix=u+1;
        for (int t=0;t<p;++t) {
            Row next{};
            int digits=0;
            long long key=0;
            for (int i=0;i<3;++i) {
                const int v=state[i]+t*row[i];
                digits+=v%p;
                next[i]=v/p;
                key=key*radix+next[i];
            }
            if (digits<p && seen.insert(key).second) todo.push_back(next);
        }
    }
    return false;
}

int main(int argc,char**argv) {
    const int bound=argc==2?std::atoi(argv[1]):2000;
    // The bound keeps every state, slope product and factorization in int range.
    if (bound<3 || bound>2000) return 2;
    for (int x=1;3*x<=bound;++x) for (int y=x;x+2*y<=bound;++y) {
        const int xy=std::gcd(x,y);
        for (int z=y;x+y+z<=bound;++z) {
            if (std::gcd(xy,z)>1) continue;
            const int u=x+y+z;
            const int H=xy*std::gcd(x,z)*std::gcd(y,z);
            if (H<=u) continue;
            const Row row{x,y,z};
            std::vector<int> ds{1};
            for (const auto &[p,w]:powers(H)) {
                int q=w;
                while (q>1 && !possible(row,p,q)) q/=p;
                const size_t old=ds.size();
                for (int v=p;v<=q;v*=p)
                    for (size_t i=0;i<old;++i) ds.push_back(v*ds[i]);
            }
            for (int t:ds) if (t>u)
                std::cout<<x<<' '<<y<<' '<<z<<' '<<t<<'\n';
        }
    }
}
