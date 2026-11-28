// Carry-free positive rows on r*x+s*y=t*z for all primitive normals in a box.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

struct Digit {
    int sum;
    bool uses_missing;
};

bool shifted_power_line(int r, int s, int t, int p) {
    if (s != r+1 || t % s != 0) return false;
    auto power = [p](int n) {
        if (n <= 1) return false;
        while (n % p == 0) n /= p;
        return n == 1;
    };
    return power(s) && power(1+t/s);
}

// One edge is already known to have a carry-free solution. Its missing
// coordinate only needs to occur once on a closed walk from zero.
bool closed_walk(int r, int s, int t, int p, bool need_x) {
    vector<vector<Digit>> by_residue(p);
    for (int x=0; x<p; ++x) for (int y=0; x+y<p; ++y) {
        for (int z=0; x+y+z<p; ++z) {
            int sum = r*x+s*y-t*z;
            bool used = need_x ? x>0 : y>0;
            if (sum==0 && used) return true;
            by_residue[(sum%p+p)%p].push_back({sum,used});
        }
    }
    int width = s+t+1;
    vector<bool> seen(2*width);
    vector<pair<int,bool>> queue{{0,false}};
    seen[t] = true;
    for (size_t i=0; i<queue.size(); ++i) {
        auto [c, used] = queue[i];
        for (Digit d : by_residue[(p-c%p)%p]) {
            int n = (c+d.sum)/p;
            bool marked = used || d.uses_missing;
            if (n==0 && marked) return true;
            int key = n+t+marked*width;
            if (!seen[key]) {
                seen[key] = true;
                queue.emplace_back(n,marked);
            }
        }
    }
    return false;
}

int main(int argc, char **argv) {
    if (argc!=3) return 2;
    int p=atoi(argv[1]), height=atoi(argv[2]);
    if (p<2 || height<1) return 2;
    for (int i=2; i<=p/i; ++i) if (p%i==0) return 2;
    long long primitive=0, weight=0, edge=0, search_free=0, family=0, other=0;
    for (int r=1; r<=height; ++r) for (int s=r; s<=height; ++s) {
        int common = gcd(r,s);
        for (int t=1; t<=height; ++t) {
            if (gcd(common,t)>1) continue;
            ++primitive;
            if ((r+t)%p==0 && (s+t)%p==0) { ++weight; continue; }
            bool x_blocked = ((r+t)/gcd(r,t))%p==0;
            bool y_blocked = ((s+t)/gcd(s,t))%p==0;
            if (!x_blocked && !y_blocked) { ++edge; continue; }
            if (x_blocked && y_blocked) return 3;
            bool has_witness = closed_walk(r,s,t,p,x_blocked);
            bool predicted = shifted_power_line(r,s,t,p);
            if (predicted && has_witness) return 4;
            if (has_witness) ++search_free;
            else if (predicted) ++family;
            else {
                if (other<30) cout << "additional " << r << ' ' << s << ' ' << t << '\n';
                ++other;
            }
        }
    }
    cout << "p=" << p << " H=" << height << " primitive=" << primitive
         << " weight=" << weight << " two_edges=" << edge << " search_free=" << search_free
         << " shifted_power=" << family << " additional=" << other << '\n';
    return 0;
}
