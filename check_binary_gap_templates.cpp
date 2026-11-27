// Mark the proportional upper and lower triples allowed by the determinant theorem.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using Triple = std::array<int,3>;

bool disjoint(const Triple& parts) {
    return !(parts[0]&parts[1]) && !(parts[0]&parts[2]) && !(parts[1]&parts[2]);
}

void mark_directions(int d, int limit, std::vector<unsigned char>& mask) {
    for (int a=(d-1)&d; a; a=(a-1)&d) {
        int rest=d^a;
        for (int b=(rest-1)&rest; b; b=(b-1)&rest) {
            int c=rest^b;
            if (!(a<b && b<c)) continue;
            int g=std::gcd(a,std::gcd(b,c));
            Triple z{a/g,b/g,c/g};
            int sum=d/g;
            for (int m=0; m*sum<=limit; ++m) {
                int low=m*sum;
                Triple parts{m*z[0],m*z[1],m*z[2]};
                if (low && disjoint(parts)) mask[low]|=1;
                if (low==limit) continue;
                for (int i=0; i<3; ++i) {
                    if (std::gcd(z[(i+1)%3],z[(i+2)%3])<3) continue;
                    ++parts[i];
                    if (disjoint(parts)) mask[low+1]|=2;
                    --parts[i];
                }
            }
        }
    }
}

int main(int argc,char** argv) {
    if (argc!=2) return 2;
    int limit=std::atoi(argv[1]);
    if (limit<1 || limit>32767) return 2;
    std::array<int64_t,4> counts{};
    for (int d=1; d<=limit; ++d) {
        std::vector<unsigned char> mask(limit+1);
        mark_directions(d,limit,mask);
        for (int c=1;c<=limit;++c) ++counts[mask[c]];
    }
    std::cout<<limit;
    for (auto total:counts) std::cout<<' '<<total;
    std::cout<<'\n';
}
