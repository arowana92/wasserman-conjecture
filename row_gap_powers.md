# Binary powers in the equality case for row gaps

Suppose a counterexample has positive gaps $d\le e$ in its two nonendpoint rows, $de=tN$ and $e-d=t-1$. For each fixed $t$ there is an effective $\delta_t>0$ such that the number of such weights $N\le X$ is
$$O_t\big((\log(2X))^{1-\delta_t}\big).$$
At $t=2661$ the equality case is excluded for every $N\le 2^{1{,}000{,}001}$. The [unrestricted-weight result](row_gap_products.md) handles every prime $t$ and every $t\le2660$.

## Two nearly equal parts in an odd base

Fix an odd prime $p$, real $\mu>0$, $K>0$, and $0\le\theta<1$. Count exponents $0\le h\le H$ for which there exist positive integers $a,b$ adding without carries in base $p$, both within $K2^{\theta h}$ of $\mu2^h$. This count is $O_{p,\mu,K,\theta}(H^{1-\delta_p})$ for some effective $\delta_p>0$ depending only on $p$.

Put $\alpha=\log 2/\log p$. The lower bound for a nonzero linear form in logarithms gives $\|k\alpha\|\gg_p k^{-C_p}$ for $k\ge1$ with an effective $C_p\ge1$. Apply it to $k\log 2-l\log p$ with $l$ the nearest integer to $k\alpha$. The form is nonzero by unique factorization, and $|l|\ll_p k$. [Baker and Wüstholz](https://doi.org/10.1515/crll.1993.442.19)

For any interval $I$ on the unit circle and any initial phase $\beta$, the number of $h\le H$ with $\{\beta+h\alpha\}\in I$ is
$$O_p\big(H|I|+H^{1-\eta}\big),\qquad \eta=1/(2C_p).$$
Indeed Dirichlet approximation at $Q=\lfloor\sqrt H\rfloor$ gives a reduced $a'/q$ with $1\le q\le Q$ and $|q\alpha-a'|\le1/Q$. The lower bound gives $q\gg_p Q^{1/C_p}$. In each block of $q$ consecutive indices the phases are within $1/Q$ of a translated grid of step $1/q$. This gives at most $q|I|+O(1)$ points, including for a partial block. There are at most $1+H/q$ blocks.

Now write $Y=\mu2^h$, $k_h=1+\lfloor\log_pY\rfloor$ and $\xi=Y/p^{k_h}\in[1/p,1)$. It equals $p^{\{\log_p\mu+h\alpha\}-1}$. An interval of length $l$ in this range pulls back to a phase interval of length $O_p(l)$ since the inverse map has bounded derivative. Put $m=\lfloor(\eta/2)\log_p H\rfloor$. For $h\ge C\log(2H)$ with a sufficiently large constant depending on the fixed parameters, both $a/p^{k_h},b/p^{k_h}$ differ from $\xi$ by less than $H^{-2}$, and $k_h\ge m$.

Within distance $H^{-2}$ of any of the $O(p^m)$ grid points $i/p^m$, count all exponents using the interval bound. Away from those neighborhoods, both normalized integers have the same first $m$ fractional digits as $\xi$. Each of these digits must be at most $(p-1)/2$, since the two integers add without a carry. This leaves at most $v^m$ intervals of length $p^{-m}$, where $v=(p+1)/2<p$. The total count, including the initial exponents, is
$$O_{p,\mu,K,\theta}\big(\log(2H)+H(v/p)^m+p^m H^{1-\eta}+p^m/H\big).$$
One can take $\delta_p=(\eta/2)(1-\log_p v)>0$.

## Applying the digit bound

Use the notation of the [row-gap proof](row_gap_products.md#when-the-right-hand-side-of-2-vanishes). At fixed $t$ only finitely many row parameters are possible apart from the two bins on the arithmetic progression. Outside the finite cases, $C$ is the row assigned 2 and $2n\mid j$. Since $w_3\mid n$ and $\gcd(j,w_3)=1$, this forces $w_3=1$. Thus the side with all three bins greater than one must be $A$. In particular $s_3>1$ is odd, because $C$ was assigned 2. Choose a prime $p\mid s_3$. The entire row $A$ adds without carries in base $p$.

Set $\lambda=c/n$, an integer, and let $P=2^h$ be the largest bit of $C_1$. Disjointness gives $C_2<P\le C_1$. Hence $C_1=P+O_t(e)$ and $N=\lambda P+O_t(e)$, with $e=O_t(\sqrt N)$ from $e=d+t-1$ and $de=tN$. The selected parts of $A$ satisfy
$$A_2=(m/b)(N-d),\qquad A_1=A_2+d.$$
Both are within $O_t(2^{h/2})$ of $(m\lambda/b)2^h$. The odd-base lemma applies with $\theta=1/2$. Each $P$ allows only $O(1)$ weights for each fixed parameter choice by the progression step in the [counting proof](row_gap_products.md#counting-the-possible-weights). Taking the minimum $\delta_p$ over the finite list proves the result for fixed $t$. When the list is empty the exceptional set is finite.

## The two remaining progressions at 2661

There are 20 parameter choices at $t=2661$. The finite checks produce twelve positive instances, of which both instances with coprime scales and bins fail the binary test. There remain two unbounded progressions, both with $(b,c,m,h',n,j)=(20,42,3,14,1,40)$ and $(s_2,s_3,w_1,w_3)=(7,3,40,1)$. Here $h'$ denotes the third multiplier in row $A$. The pairs $(x,y)$ are $(3,887)$ and $(2661,1)$.

In general on such a progression set $K=t\lambda$. Its binary row is determined by $e$ alone:
$$C_1=e(e+1)/K,\qquad C_2=C_1-e,\qquad C_3=(\lambda-2)C_1.$$
Indeed $\lambda C_1=N+e=e(e+1)/t$. The inequalities $C_2<P\le C_1$ imply
$$e(e+1-K)<KP\le e(e+1).$$
If $s=\lfloor\sqrt{KP}\rfloor$, then $s\le e\le s+K$. The lower bound follows from $(s-1)s<KP$. For the upper bound, $e\ge s+K+1$ would give $e(e+1-K)>(s+1)^2>KP$.

Here $K=111762$. The congruences $e\equiv-1\bmod xc$ and $e\equiv t\bmod yb$ give the two residues $357461$ and $111761$, respectively, modulo $L=1117620$. For $P=2^h$ and $0\le h\le1{,}000{,}000$, the intervals $[s,s+K]$ give 99,989 and 99,656 progression points, respectively. Every one produces a binary carry among the first 40 digits of $C$. Computing $C_1\bmod 2^{40}$ needs only $e\bmod K2^{40}$ by the exact division formula.

One program obtains the bits of $\lfloor\sqrt{K2^h}\rfloor$ from a single integer square root for each parity of $h$. Shifting its high-precision root by one bit gives the root at the preceding exponent of that parity. Independently, another program computes $z=\lfloor\sqrt{K2^h}/L\rfloor$ and tests both $e=Lz+r$ and $e=L(z+1)+r$ for each possible residue $r$. These larger windows contain every candidate since $K<L$. This second calculation checked 4,000,004 points. All failed its binary digit-sum test. Both programs square their top roots to verify the exact floors. Since $C_1<N$, every weight $N\le2^{1{,}000{,}001}$ falls within the exponent bound.

Run the following with Python 3.9 or later. The last command checks the first 4,097 exponents with a fresh integer square root at each one.

~~~sh
python3 check_gap_dyadic_windows.py
python3 check_gap_dyadic_quotients.py
python3 check_gap_dyadic_windows.py 4096 --direct
~~~
