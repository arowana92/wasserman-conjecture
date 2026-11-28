// List multipliers meeting the two-base tests using factorial-floor intervals.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using U = std::uint64_t;

int main(int argc, char** argv) {
    if (argc!=5) return 2;
    U a=std::strtoull(argv[1],nullptr,10), b=std::strtoull(argv[2],nullptr,10);
    U limit=std::strtoull(argv[3],nullptr,10), first=std::strtoull(argv[4],nullptr,10);
    if (!((a>=2 && a<=4 && b==5) || (a==5 && b==6)) ||
        (first!=2 && first!=3) || limit>1000000000000000000ULL/(a+b)) return 2;
    std::vector<U> layers;
    U modulus=1;
    for (U p : {first,U(5)}) {
        U w=(a%p==0?a:b), power=1;
        if (w%p!=0) return 2;
        while (w%p==0) { w/=p; power*=p; }
        modulus*=power;
        for (U q=p; q<=(a+b)*limit;) {
            layers.push_back(q);
            if (q>(a+b)*limit/p) break;
            q*=p;
        }
    }
    std::sort(layers.rbegin(),layers.rend());
    U residue=0;
    for (U r=0;r<modulus;++r) {
        bool ok=true;
        for (U p : {first,U(5)}) {
            U v=(a%p==0?a:b), unit=(a%p==0?b:a), q=1;
            while (v%p==0) { v/=p; q*=p; }
            if ((unit*r)%q!=1) ok=false;
        }
        if (ok) residue=r;
    }
    if (!residue) return 3;
    for (U h=residue; h<=limit;) {
        U x=a*h, y=b*h, next=h;
        for (U q : layers) {
            U rx=x%q, ry=y%q;
            if (rx+ry<q) continue;
            // The carry persists until one of the two part quotients increases.
            U dx=(q-rx+a-1)/a, dy=(q-ry+b-1)/b;
            next=h+std::min(dx,dy);
            break;
        }
        if (next==h) {
            std::cout << h << '\n';
            h+=modulus;
        } else {
            h=next+(residue+modulus-next%modulus)%modulus;
        }
    }
}
