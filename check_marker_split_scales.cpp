// Recover every candidate small row scale at the largest binary marker.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using I = int64_t;
using Wide = __int128_t;

I inverse(I a, I modulus) {
    I r = a, s = modulus, x = 1, y = 0;
    while (s) {
        I k = r/s, next_r = r-k*s, next_x = x-k*y;
        r = s; s = next_r; x = y; y = next_x;
    }
    return r == 1 ? (x%modulus+modulus)%modulus : 0;
}

bool binary_split(I u,I t) {
    std::vector<I> possible;
    for (I x=u; x<t; x+=u) {
        if ((x&t)!=x) continue;
        for (I y:possible)
            if (!(x&y) && x+y<t) return true;
        possible.push_back(x);
    }
    return false;
}

int main() {
    constexpr I r=(I(1)<<31)-1, high=r+1, dmax=512, jmax=2, kmax=6;
    constexpr I amax=dmax<<jmax;
    I scales=0,divisors=0;
    for (I z=1; z<amax; ++z) {
        I iz=inverse(z,r);
        if (!iz) return 1;
        for (I a=1; 4*a<dmax; ++a) {
            if (std::gcd(a,z)!=1) continue;
            I residue=(a*iz)%r;
            for (I block=2*a; block<dmax/2; ++block) {
                I u=block*r+residue;
                if (u%2==0) continue;
                ++scales;
                I hi_inverse=inverse(high,u);
                I m=u-inverse(r,u);
                for (I k=1; k<=kmax; ++k) {
                    m=I(Wide(m)*hi_inverse%u);
                    if (m>amax) continue;
                    for (I j=0; j<=jmax; ++j) {
                        I shift=I(1)<<j, d=m/shift;
                        if (d*shift!=m || !d || d>dmax || 2*u>=d*r) continue;
                        ++divisors;
                        if (binary_split(u,d*r))
                            std::cout<<d<<' '<<k<<' '<<j<<' '<<u<<'\n';
                    }
                }
            }
        }
    }
    std::cerr<<scales<<" candidate scales, "<<divisors<<" scale and weight incidences\n";
}
