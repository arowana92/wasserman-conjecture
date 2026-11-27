// A necessary binomial-cut test for the row with predecessor product >= its normalized sum.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using U64 = std::uint64_t;
using U128 = __uint128_t;
struct Counts { U64 weights = 0, subsets = 0, cuts = 0, rows = 0; };
struct Factor { U128 prime, power; };
// The allowance is the number of valuation units removed from a full power.
struct SmallPower { unsigned prime, power; int allowance=0; };
constexpr unsigned bound = 67;
constexpr unsigned largest_input_power = 1000000;
unsigned inverse[bound+1][bound+1];

unsigned inverse_mod(unsigned a,unsigned w) {
    if(w<=bound) return inverse[w][a];
    std::int64_t r=w,s=a,t=0,z=1;
    while(s) {
        auto k=r/s, next=r-k*s, tz=t-k*z;
        r=s; s=next; t=z; z=tz;
    }
    if(r!=1) std::abort();
    return unsigned((t%w+w)%w);
}

std::string decimal(U128 n) {
    std::string s;
    do { s.push_back(char('0'+n%10)); n/=10; } while(n);
    std::reverse(s.begin(),s.end());
    return s;
}
U128 number(const std::string& s) {
    if(s.empty()) std::abort();
    U128 n=0;
    for(char c:s) {
        if(c<'0' || c>'9' || n>(~U128(0)-unsigned(c-'0'))/10) std::abort();
        n=10*n+unsigned(c-'0');
    }
    return n;
}
template<typename Z>
bool cut_without_carry(Z n,Z a,unsigned p) {
    if(p==2) return (a & (n-a))==0;
    while(a) {
        if(a%p>n%p) return false;
        a/=p; n/=p;
    }
    return true;
}

template<typename Z>
int cut_valuation(Z n,Z part,unsigned p,int limit) {
    if(limit==0) return cut_without_carry(n,part,p) ? 0 : 1;
    int total=0;
    Z other=n-part;
    while(n && total<=limit) {
        part/=p; other/=p; n/=p;
        total+=int(n-part-other);
    }
    return total;
}

template<typename Z>
struct Search {
    Z n;
    const std::vector<SmallPower>& factors;
    bool retained;
    Counts& count;
    std::vector<SmallPower> selected;
    std::vector<Z> residues;
    std::vector<int> admissible;

    void cuts(Z u,Z v,Z modulus) {
        ++count.subsets;
        residues.clear();
        residues.push_back(0);
        // Each subset chooses the exceptional residues of one normalized part.
        for(auto [p,w,allowance]:selected) {
            Z m=modulus/w;
            unsigned multiplier=unsigned(U64(v%w)*inverse_mod(unsigned(m%w),w)%w);
            Z unit=m*multiplier;
            int size=int(residues.size());
            for(int i=0;i<size;++i) {
                Z x=residues[i];
                residues.push_back(x>=modulus-unit ? x-(modulus-unit) : x+unit);
            }
        }
        admissible.clear();
        int full=int(residues.size())-1;
        // A smallest part has size at most one third of the normalized sum.
        for(int mask=1;mask<full;++mask) {
            Z x=residues[mask];
            if(!x || x>v/3) continue;
            Z part=u*x;
            bool ok=true;
            for(auto [p,w,allowance]:selected) {
                // The discarded digits are 0 or 1 in the part and 1 in n.
                if(cut_valuation<Z>(n/w,part/w,p,allowance)>allowance) { ok=false; break; }
            }
            if(ok) admissible.push_back(mask);
        }
        count.cuts+=admissible.size();
        for(int mask:admissible) {
            Z x=residues[mask];
            int rest=full^mask;
            for(int j=(rest-1)&rest;j;j=(j-1)&rest) {
                int k=rest^j;
                if(j>=k) continue;
                Z y=residues[j], z=residues[k];
                if(y<x || z<x || y>=v-x || z!=v-x-y) continue;
                bool ok=true;
                // Factor the whole coefficient into the cut and the other two parts.
                for(auto [p,w,allowance]:selected) {
                    int first=cut_valuation<Z>(n/w,(u*x)/w,p,allowance);
                    if(cut_valuation<Z>(n-u*x,u*y,p,allowance-first)>allowance-first) {
                        ok=false; break;
                    }
                }
                if(!ok) continue;
                ++count.rows;
                std::cout<<decimal(n)<<' '<<decimal(u)<<' '<<decimal(modulus)
                         <<' '<<decimal(u*x)<<' '<<decimal(u*y)<<' '<<decimal(u*z)<<'\n';
            }
        }
    }

    void visit(Z u,Z v,Z product,int start) {
        // Stop exactly when the selected prefix determines all parts uniquely.
        if(product>=v) {
            if(selected.size()>=3) cuts(u,v,product);
            return;
        }
        for(int i=start;i<int(factors.size());++i) {
            auto power=factors[i];
            while(power.power>1) {
                if(power.power<=(v-1)/2) {
                    selected.push_back(power);
                    visit(u,v,product*power.power,i+1);
                    selected.pop_back();
                }
                if(!retained) break;
                power.power/=power.prime;
                ++power.allowance;
            }
        }
    }
    void run(const std::vector<Factor>& weight) {
        std::vector<Z> scales{1};
        for(const auto& f:weight) {
            int len=int(scales.size());
            Z w=Z(f.power);
            while(w>1) {
                for(int i=0;i<len;++i) scales.push_back(scales[i]*w);
                if(!retained) break;
                w/=Z(f.prime);
            }
        }
        for(Z u:scales) {
            if(u >= (retained ? 2u : 3u) && n/u >= 3) visit(u, n/u, 1, 0);
        }
    }
};

int main(int argc,char** argv) {
    bool retained=argc==2 && std::string(argv[1])=="--retained";
    if(argc!=1 && !retained) return 2;
    for(unsigned w=2;w<=bound;++w)
        for(unsigned i=1;i<w;++i)
            for(unsigned j=1;j<w;++j) if(i*j%w==1) inverse[w][i]=j;
    Counts counts;
    std::string word;
    while(std::cin>>word) {
        U128 n=number(word);
        int k;
        if(n<3 || !(std::cin>>k) || k<1 || k>40) return 2;
        std::vector<Factor> nf(k);
        U128 product=1;
        for(auto& f:nf) {
            std::string p,w;
            if(!(std::cin>>p>>w)) return 2;
            f={number(p),number(w)};
            if(f.prime<2 || f.power<f.prime || product>n/f.power) return 2;
            product*=f.power;
        }
        if(product!=n || !(std::cin>>k) || k<1 || k>20) return 2;
        std::vector<SmallPower> mf(k);
        product=1;
        for(auto& f:mf) {
            if(!(std::cin>>f.prime>>f.power) || f.power>largest_input_power || f.prime<2 || f.power<f.prime || product>(n-1)/f.power) return 2;
            product*=f.power;
        }
        if(product!=n-1) return 2;
        std::sort(mf.begin(),mf.end(),[](auto a,auto b){return a.power>b.power;});
        ++counts.weights;
        if(n<=~U64(0)) {
            Search<U64> s{U64(n),mf,retained,counts,{},{},{}};
            s.run(nf);
        } else {
            Search<U128> s{n,mf,retained,counts,{},{},{}};
            s.run(nf);
        }
    }
    std::cerr<<counts.weights<<" weights, "<<counts.subsets<<" subsets, "
             <<counts.cuts<<" smallest cuts, "<<counts.rows<<" rows\n";
}
