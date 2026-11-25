# A quantitative common divisor from two parts of one row

Write $K_A$ for the trinomial coefficient of a positive row $A$ of weight $N$. Let $A,C$ be **any** two such rows and let $a,b$ be **any** two parts of $A$. Put $s=a+b$, $D=\gcd(N,K_A,K_C)$ and $J=\gcd(N-1,K_A,K_C)$. Then
$$729(N-1)^4\le D^6J^4(ab)^3(s-1).\tag{1}$$
In particular the common divisor of $K_A,K_C$ and the endpoint coefficient $N(N-1)$ satisfies
$$DJ\ge\frac{3(N-1)^{2/3}}{\sqrt{ab}\,(s-1)^{1/6}}.$$
For fixed $a,b$ this lower bound tends to infinity uniformly in the other row $C$.

**Endpoint consequence.** In a counterexample containing $E=(1,1,N-2)$, every pair of parts of either other row satisfies
$$s^7>6^6(N-1)^4.\tag{2}$$
Indeed $D=J=1$, $ab\le s^2/4$, and $s-1<s$. In particular the exponent for a sum of two parts is at least $4/7$, even though the earlier [individual-part bound](endpoint_factors.md) has exponent $2/5$.

## Retaining partial prime powers

For $p^e\parallel N$ or $p^e\parallel N-1$, let $h=0$ or $1$, respectively. For each row the valuation of its coefficient is the sum over $j\ge1$ of
$$\frac{\sum_i (A_i\bmod p^j)-(N\bmod p^j)}{p^j}.$$
This follows by subtracting Legendre’s three factorial valuations. The summands are nonnegative integers. If the low residues at one level $j\le e$ sum to more than $h$, their sum also exceeds $h$ at every level from $j$ through $e$, since residues can only increase and $N\bmod p^j=h$ throughout. Thus if $v_p(K_A)=l<e$, the residues at level $e-l$ sum to $h$.

Retain $p^{e-\min(e,v_p(K_A),v_p(K_C))}$ and assign it to a row with smaller valuation, choosing either on a tie. The products assigned to $A,C$ from $N$ are $U,V$, and those from $N-1$ are $S,T$. The removed products are exactly $D,J$, so
$$DUV=N,\qquad JST=N-1.$$
The parts of $A$ are multiples of $U$, and those of $C$ are multiples of $V$. At each retained power from $N-1$ the three residues of its assigned row are $1,0,0$ in some order.

Let $x=a/U$, $y=b/U$. At a retained power $w$ in $S$ at least one of $x,y$ is a multiple of $w$. If both are multiples, $w^2\mid xy$. Otherwise the actual residues of $a,b$ sum to 1, giving $w\mid s-1$. Consequently
$$S\mid xy,\qquad S^2\mid xy(s-1),\qquad S^4\mid (xy)^3(s-1).$$
For the other row write its normalized parts as $z_i=C_i/V$. At each retained power in $T$, two $z_i$ are multiples of that power since $V$ is coprime to $T$. Hence
$$T^2\mid z_1z_2z_3\le (DU/3)^3.$$
Combining these inequalities cancels $U^6$ and proves (1).

Two separate consequences when $D=J=1$ will also be useful. Write $g=\gcd(a,b)$. Since $S\mid\operatorname{lcm}(x,y)=ab/(Ug)$, and $U\mid g$,
$$ab\ge 3\sqrt3\,(N-1)g/\sqrt U\ge3\sqrt3\,(N-1)\sqrt g.\tag{3}$$
Using $S^2\mid xy(s-1)$ instead gives
$$27(N-1)^2\le ab(s-1)U.\tag{4}$$

## A third row with a specified part

Suppose $K_A,K_B,K_C$ have gcd one. Every prime in $DJ$ must be avoided by $B$. If $m\ge2$ is a part of $B$, this gives $D\mid m$ and $J\mid m(m-1)$ by the low-residue conditions. Since $D,J$ are coprime, $J\mid(m/D)(m-1)$. In particular
$$s^7 m^6(m-1)^4>6^6(N-1)^4.\tag{5}$$
This applies to every pair of parts in a different row, for arbitrary $N$.

If $B=(1,v,w)$, then $D=1$ and $J\mid\gcd(v,w)$. In this case
$$\gcd(v,w)>6^{3/2}(N-1)/s^{7/4}.\tag{6}$$
For example, take $A$ to be a row avoiding the largest prime at most $N$. Its two smaller parts sum to $O(N^{21/40})$ by Baker–Harman–Pintz. Thus a part $m\ge2$ in **any other row** must satisfy $m\gg N^{13/400}$, and a 1 in such a row requires its two companions to have gcd $\gg N^{13/160}$. [Baker–Harman–Pintz](https://doi.org/10.1112/plms/83.3.532)

## Reproduction

Run `python3 check_endpoint_pair_bound.py`. It forms exact binomial products for every positive unordered row of weights 3 through 45 and checks every ordered choice of two rows and every choice of two positions in the first. It tests the retained divisibilities as well as (1). Among 810,876 checks, 123,738 include a partially retained prime power and 795,030 yield a lower bound strictly greater than one for the common divisor.
