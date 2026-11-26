// Independently reconstruct V first, then check individual entries and factorial floors.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using Z=int64_t;
using Row=std::array<Z,3>;
using Powers=std::vector<std::pair<int,int>>;

int max_t(int u) {
    int low=0, high=200000000;
    const __int128 cube=__int128(u)*u*u;
    while (low+1<high) {
        int m=low+(high-low)/2;
        if (__int128(27)*m*m<=cube) low=m;
        else high=m;
    }
    return low;
}

struct FactorTable {
    std::vector<uint16_t> factor;
    explicit FactorTable(int n) : factor(n/2+1) {
        for (int p=3; Z(p)*p<=n; p+=2) {
            if (factor[p/2]) continue;
            for (int j=p*p; j<=n; j+=2*p) factor[j/2]=p;
        }
    }
    Powers powers(int n) const {
        Powers list;
        while (n>1) {
            int p=n%2 ? factor[n/2] : 2;
            if (!p) p=n;
            int w=1;
            do { n/=p; w*=p; } while (n%p==0);
            list.emplace_back(p,w);
        }
        return list;
    }
};

std::vector<Z> all_divisors(const Powers& factors) {
    std::vector<Z> out{1};
    for (auto [p,w] : factors) {
        auto count=out.size();
        for (int z=w; z>1; z/=p)
            for (size_t i=0;i<count;++i) out.push_back(out[i]*z);
    }
    return out;
}

// Extended Euclid by backward substitution, with coprime inputs below 200 million.
Z inverse_fast(Z a, Z m) {
    if (m==1) return 0;
    if (a==1) return 1;
    if (!a) std::abort();
    Z t=inverse_fast(m%a,a);
    return m-(m*t-1)/a;
}

bool binomial_unit(Z n, Z a, int p) {
    if (p==2) return __builtin_popcountll(a)+__builtin_popcountll(n-a)==__builtin_popcountll(n);
    Z b=n-a;
    while (a) {
        n/=p; a/=p; b/=p;
        if (n!=a+b) return false;
    }
    return true;
}

bool trinomial_unit(Z n, Row row, int p) {
    do {
        n/=p;
        for (Z& a : row) a/=p;
        if (n!=row[0]+row[1]+row[2]) return false;
    } while (n);
    return true;
}

void print_row(char side, int a, int u, Z s, int t, Z n, Row row) {
    std::sort(row.begin(),row.end());
    std::cout << side << ' ' << a << ' ' << u << ' ' << s << ' ' << t << ' ' << n;
    for (Z b : row) std::cout << ' ' << b;
    std::cout << '\n';
}

int audit_companion(int a,int u,Z s,int t,Z n,const Powers& pt) {
    // For a unitary divisor h of T, a coordinate is 0 at T/h and U at h.
    std::vector<int> unions{1};
    for (auto [p,w] : pt) {
        size_t old=unions.size();
        for (size_t j=0;j<old;++j) unions.push_back(unions[j]*w);
    }
    std::vector<int> entries;
    for (int h : unions) {
        if (h==1) continue;
        int q=t/h;
        if (q>=u) continue;
        Z y=q*(Z(u%h)*inverse_fast(q%h,h)%h);
        if (0<y && y<u) entries.push_back(int(y));
    }
    std::sort(entries.begin(),entries.end());
    int found=0;
    for (size_t i=0; i<entries.size(); ++i) {
        for (size_t j=i; j<entries.size(); ++j) {
            int z=u-entries[i]-entries[j];
            if (z<entries[j]) break;
            if (!std::binary_search(entries.begin()+j,entries.end(),z)) continue;
            Z v=n/u;
            Row row{v*entries[i],v*entries[j],v*z};
            bool good=true;
            for (auto [p,w] : pt) if (!trinomial_unit(n,row,p)) good=false;
            if (good) { print_row('C',a,u,s,t,n,row); ++found; }
        }
    }
    return found;
}

int audit_anchored(int a,int u,Z s,int t,Z n,int b,const Powers& pu,const Powers& ps) {
    int found=0;
    Z step=Z(b)*u;
    for (Z z=step; 2*z<=n-a; z+=step) {
        Row row{a,z,n-a-z};
        bool good=true;
        for (auto [p,w] : ps) if (!trinomial_unit(n,row,p)) { good=false; break; }
        if (!good) continue;
        for (auto [p,w] : pu) if (!trinomial_unit(n,row,p)) { good=false; break; }
        if (good) { print_row('A',a,u,s,t,n,row); ++found; }
    }
    return found;
}

void audit(int first,int last) {
    FactorTable tab(std::max(last,max_t(last)));
    Z systems=0,heavy_a=0,heavy_c=0,rows_a=0,rows_c=0,largest=0;
    for (int u=3;u<=last;++u) {
        const auto pu=tab.powers(u);
        int hi=max_t(u);
        for (int x=(first+u-1)/u; x<=last/u; ++x) {
            int a=u*x;
            auto factors=tab.powers(x), pred=tab.powers(a-1);
            factors.insert(factors.end(),pred.begin(),pred.end());
            for (Z s : all_divisors(factors)) {
                if (std::gcd(s,Z(u))!=1) continue;
                int b=int(std::gcd(s,Z(a-1)));
                Z min_v=x+2LL*b;
                Z v0=(1+s*(u-inverse_fast(s%u,u)))/u;
                Z v=v0;
                if (v<min_v) v+=((min_v-v+s-1)/s)*s;
                if (u*v>s*hi+1) continue;
                Powers ps;
                for (auto [p,w] : factors) if (s%p==0) ps.emplace_back(p,0);
                for (; u*v<=s*hi+1; v+=s) {
                    Z n=u*v;
                    int t=int((n-1)/s);
                    if (std::gcd(v,Z(u))!=1 || std::gcd(s,Z(t))!=1) continue;
                    ++systems; largest=std::max(largest,n);
                    bool ok=true;
                    for (auto [p,w] : pu) if (!binomial_unit(n,a,p)) { ok=false; break; }
                    if (!ok) continue;
                    for (auto [p,w] : ps) if (!binomial_unit(n,a,p)) { ok=false; break; }
                    if (!ok) continue;
                    auto pt=tab.powers(t);
                    for (auto [p,w] : pt) if (w>=(u+1)/2) ok=false;
                    if (!ok) continue;
                    if (s>v) {
                        ++heavy_a;
                        rows_a+=audit_anchored(a,u,s,t,n,b,pu,ps);
                    } else if (pt.size()>=3) {
                        ++heavy_c;
                        rows_c+=audit_companion(a,u,s,t,n,pt);
                    }
                }
            }
        }
    }
    std::cerr << "range " << first << ' ' << last << " systems " << systems << " largest "
              << largest << " heavy A/C " << heavy_a << ' ' << heavy_c
              << " passing rows A/C " << rows_a << ' ' << rows_c << '\n';
}

int main(int argc,char** argv) {
    if (argc==2 && std::string(argv[1])=="--example") {
        FactorTable tab(887603);
        const int n=222788354,s=887603,u=283;
        for (int a : {86617810,45491401,90679143})
            audit_anchored(a,u,s,251,n,int(s/std::gcd(a/u,s)),tab.powers(u),tab.powers(s));
        audit_companion(787238,787238,251,s,n,tab.powers(s));
        return 0;
    }
    if (argc!=3) return 2;
    int first=std::stoi(argv[1]),last=std::stoi(argv[2]);
    if (first<3 || first>last || last>1000000) return 2;
    audit(first,last);
}
