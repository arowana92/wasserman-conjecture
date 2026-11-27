// Project the first and second parts by enumerating the leading total.
// Keep both endpoints as exact fractions throughout the sweep.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using I = std::int64_t;
struct Fraction {
    I n,d;
    bool operator<(Fraction x) const { return n*x.d < x.n*d; }
};
struct Interval { Fraction a,b; };

void cuts(std::vector<Interval>& out, int n, int difference, int rest,
          int b, int place, int p) {
    if (rest == 0) {
        out.push_back({{b,n+1},{5*b+2+difference,5*(n+1)}});
        return;
    }
    for (int j = 0; j <= rest%p; ++j)
        cuts(out,n,difference,rest/p,b+j*place,place*p,p);
}

std::vector<Interval> cover(int p, int height) {
    int limit = 1;
    for (int i=0; i<height; ++i) limit *= p;
    std::vector<Interval> raw;
    for (int n=limit/p; n<limit; ++n) {
        int a = (3*n+2)/5;
        bool ok = true;
        for (int x=n,y=a; x; x/=p,y/=p) {
            int difference = 5*y-3*x;
            if (difference < -2 || difference > 2 || y%p > x%p) {
                ok = false;
                break;
            }
        }
        if (ok) cuts(raw,n,5*a-3*n,n-a,0,1,p);
    }
    std::sort(raw.begin(),raw.end(),[](Interval x, Interval y){return x.a<y.a;});
    std::vector<Interval> result;
    for (Interval x : raw) {
        if (result.empty() || result.back().b<x.a) result.push_back(x);
        else result.back().b = std::max(result.back().b,x.b);
    }
    return result;
}

std::vector<Interval> intersect(const std::vector<Interval>& a,
                                const std::vector<Interval>& b) {
    std::vector<Interval> result;
    std::size_t i=0, j=0;
    while (i<a.size() && j<b.size()) {
        Interval x{std::max(a[i].a,b[j].a),std::min(a[i].b,b[j].b)};
        if (!(x.b<x.a)) result.push_back(x);
        if (a[i].b<b[j].b) ++i;
        else ++j;
    }
    return result;
}

int main() {
    auto result = intersect(intersect(cover(2,20),cover(3,12)),cover(5,8));
    int positive = 0;
    I upper = 0;
    Fraction rightmost{0,1};
    const I grid = 1'000'000'000'000;
    for (Interval i : result) {
        if (!(i.a<i.b)) continue;
        ++positive;
        upper += (i.b.n*grid+i.b.d-1)/i.b.d - i.a.n*grid/i.a.d;
        if (i.a<Fraction{1,5}) rightmost = std::max(rightmost,i.b);
    }
    assert(result.size()==12978 && positive==11742);
    assert(upper < 15'382'000'000);
    assert(rightmost.n*2'286'935 == rightmost.d*95'574);
    std::cout << result.size() << " components, " << positive
              << " positive lengths, total length at most " << upper << '/' << grid << '\n';
}
