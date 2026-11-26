// Recover bin triples from the integer (d+e+f)-9(T-1)/U, independently of resultants.
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

using Z=int64_t;
using Wide=__int128_t;

void show(Wide n) {
    if (n>=10) show(n/10);
    std::cout << int(n%10);
}

bool square(Z n,Z& r) {
    if (n<0) return false;
    // A square can end in just these four residues modulo 16.
    if (!((0x0213 >> (n%16)) & 1)) return false;
    r=Z(std::sqrt(static_cast<long double>(n)));
    while (r*r>n) --r;
    while ((r+1)*(r+1)<=n) ++r;
    return r*r==n;
}

bool test(Z d,Z h,Z k) {
    Z e=d+h,f=d+k;
    if (d<=1 || d%2==0 || std::gcd(d,e)!=1 || std::gcd(e,f)!=1 || std::gcd(d,f)!=1)
        return false;
    Z u=d*e+d*f+e*f;
    Wide t=Wide(d)*e*f;
    if ((t-1)%u || (t-1)/u<2 || std::gcd(Z(t%(u-1)),u-1)!=1) return false;
    std::cout << "BINS " << d << ' ' << e << ' ' << f << ' ' << u << ' ';
    show((u-1)*t+1);
    std::cout << '\n';
    return true;
}

int main(int argc,char** argv) {
    if (argc!=2) return 2;
    int gap=std::stoi(argv[1]);
    if (gap<4 || gap>4096) return 2;
    Z pairs=0,triples=0;
    for (Z k=4;k<=gap;k+=2) {
        for (Z h=2;h<k;h+=2) {
            ++pairs;
            Z s=h+k,b=h*k,a=2*(h*h-h*k+k*k);
            for (Z m=3;m<=s+1;m+=4) {
                Z linear=2*m*s-a;
                Z constant=(m-s)*b-9;
                Z root;
                if (!square(linear*linear-12*m*constant,root)) continue;
                for (int sign : {-1,1}) {
                    if (sign>0 && !root) continue;
                    Z numerator=-linear+sign*root;
                    if (numerator<=0 || numerator%(6*m)) continue;
                    triples+=test(numerator/(6*m),h,k);
                }
            }
        }
    }
    std::cerr << "gap pairs " << pairs << " bin triples " << triples << '\n';
}
