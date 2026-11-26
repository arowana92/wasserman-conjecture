# Adjacent cuts near an edge

Suppose $1\le a<a+1\le N/2$, $M=N-1$ and $\gcd(NM,\binom Na,\binom N{a+1})=1$. Then
$$108a^5>M^4,\qquad 470^4N^4a^{21}>M^{21}. \tag{1}$$
In particular, the smaller index satisfies $a\gg N^{17/21}$. The first inequality and the [adjacent-pair check through one billion](adjacent_cuts.md) rule out every $a\le6{,}213{,}198$ at **every weight**.

Here is a consequence for whole rows of Pascal's triangle. Let $P=p^k\le N<2P$ with $p\mid N(N-1)$. Every adjacent pair of interior binomial coefficients of weight $N$ shares a prime factor of $N(N-1)$ whenever either
$$108(N-P)^5\le M^4,\qquad 470^4N^4(N-P)^{21}\le M^{21}.$$
Indeed, if a binomial avoids $p$, one of the two sides of its cut contains the leading base-$p$ unit $P$. Thus $N-a\ge P$ for a pair with gcd one, contradicting (1). Taking $P=2^k$ gives intervals of length $\gg 2^{17k/21}$ above every sufficiently large power of two. The result applies to the trinomial problem when the two nonendpoint rows have adjacent cuts.

## Two positive integers from the small cut

Use the [adjacent-cut parameters](adjacent_cuts.md). They satisfy
$$M=gi,\quad a=gc,\quad c+f=i,\quad cf=Nj,\quad \ell i=1+g^2j,\quad g>4\ell.$$
All parameters are positive integers, $c\le f$, and $\gcd(g,\ell)=1$. Write $x=g/\ell>4$ and $v=j\ell^2$, so $N=1+x+vx^3$. Since $(a+1)(N-a)=Ni\ell$ and $N-a<N$, integrality gives $a\ge i\ell=M/x$.

Put $k=c-gj$. Then $k=(c^2+j)/i>0$ and
$$k(f-gj)=j\ell i.$$
Since $gj=(cf-j)/i<f/2$, we have $f-gj>i/4$ and $k<4j\ell$. Rearrangement of the last product gives the positive integer $b=k-j\ell$ with
$$bi=k(2gj+k),\qquad 0<b<12j\ell^2/g.$$
The bound uses $g>4\ell$ and $i>g^2j/\ell$. Set $w=bg-\ell k$. Completing the square in $b(1+g^2j)=\ell k(2gj+k)$ gives $jw^2=\ell k^3-b^2$. For $B=2\ell^2j+3\ell b$ this reads
$$4w^2-B^2=3\ell^2b^2+4b^2(\ell b-1)/j>0.$$
Since $w>-\ell k$ and $B>2\ell k$, we get $2w>B$. Also
$$(B+3b^2/(4j))^2-4w^2=\ell b^3/(2j)+4b^2/j+9b^4/(16j^2)>0.$$
It follows that $z=2w-B$ is a positive integer satisfying
$$0<z<\frac{3b^2}{4j}<\frac{108v}{x^2}<108M/x^5.$$
Using $a\ge M/x$ proves the first inequality of (1). The other positive integer $T=3b^2-4jz$ satisfies
$$T(2\ell^2j+3\ell b)=\ell b^3+2jz^2+8b^2. \tag{2}$$
This follows by substituting $2w=B+z$ above. Since $b/(j\ell)<12/x<3$, the right side of (2) is less than $13\ell b^3$ by $z<3b^2/(4j)$. Thus
$$1\le T<7b^3/(\ell j)<12100\, v^2/(\ell^2 x^3). \tag{3}$$

## A rational approximation and its zero case

Let $y=b/(j\ell)$, $p=6x-13$ and $q=3x^2-14x+14$. Dividing $jw^2=\ell k^3-b^2$ by $j^3\ell^4$ gives
$$F(y)=0,\qquad F(Y)=Y^2-(x^2-2x-2+1/v)Y+2x+1.$$
Indeed $[(x-1)y-1]^2=(1+y)^3-y^2/v$ and $y>0$. The rational approximation $p/q$ satisfies
$$F(p/q)=(12x+1)/q^2-p/(vq).$$
For $x\ge16$, $2x^2\le q\le3x^2$, $0<p/q\le3/x$, and $y<12/x$. Hence $|F(p/q)|\le4/x^3+3/(vx)$ and $|y+p/q-(x^2-2x-2+1/v)|\ge x^2/2$. Dividing the difference $F(p/q)-F(y)$ by this last factor yields
$$|q y-p|\le 24/x^3+18/(v x).$$
The quantity $L=\ell v(qy-p)$ is an integer, explicitly
$$L=(3g^2-14g\ell+14\ell^2)b-v(6g-13\ell).$$
Using $v>x^2/108$ from the bound for $z$ gives $|L|<2000\ell v/x^3$.

If $L=0$ then $F(p/q)=0$, giving
$$j\ell^4(12g+\ell)=(6g-13\ell)(3g^2-14g\ell+14\ell^2).$$
Coprimality gives $\ell\mid18$ by reduction modulo $\ell$. Multiplying the right side by $12^3$ and reducing modulo $12g+\ell$ gives $-162\cdot2187\,\ell^3$. Enumerating the divisors of $354294\ell^3$ for the six possibilities for $\ell$ yields only $(g,\ell,j)=(20,3,2),(182,3,546)$ with $g>4\ell$ and $j$ integral and positive. The first has $x<16$ and the second fails $\ell\mid1+g^2j$.

Thus, for $x\ge16$, $v\ell>x^3/2000$. By (3), $v/\ell>x^{3/2}/110$. Multiplication gives $v>x^{9/4}/470$ and $N>x^{21/4}/470$. For $x<16$ the same last inequality holds for $N>10^9$, and the finite adjacent-cut theorem covers smaller $N$. Combining with $a\ge M/x$ proves the second inequality of (1).

Run `python3 check_adjacent_cuts.py` to reproduce the zero-case enumeration and to check the two positive integers, identity (2) and the bounds on the reconstructed pairs. The calculation of the first all-weight cutoff uses $108\cdot6213198^5\le(10^9)^4<108\cdot6213199^5$.
