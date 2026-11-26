# One four-prime list for $2^{76}$ weights

Set
$$P=2^{256},\quad Q=P+297,\quad R=P+3072559256396398118163687,\quad F=PQR,\quad r=F-1185.$$
The **same four primes $\{2,Q,R,r\}$** certify every selection of three positive trinomial coefficients at all the weights
$$dF,\quad dF+1,\qquad 1\le d\le 2^{75}.$$
These are $2^{76}=75{,}557{,}863{,}725{,}914{,}323{,}419{,}136$ distinct weights. The same list has minimum possible size at both weights for at least $16\lfloor 2^{75}/225\rfloor$ of the multipliers. In particular the number four is necessary for any one list covering the entire range. The argument proves the three-row conclusion directly by pairwise incompatibility.

## A fourth prime usable at every multiplier

Let $P=p^e$, $Q=P+a$, $R=P+b$ be powers of distinct primes with $0<2a<b$, and put $F=PQR$. Let $B=p^l$ be a power of $p$ with $QR\equiv-1\bmod B$. Suppose $r=F-g$ is a prime different from the three base primes, $g>0$, and a positive integer $D$ satisfies
$$D\le B,\qquad 2Db^2\le P,\qquad Dg+1<P,\qquad D<r.$$
Then the three bases together with $r$ give a single certificate at all $dF,dF+1$, $1\le d\le D$.

The bases of $P,Q,R$ are pairwise incompatible for trinomials at both weights by the [three-power result](clustered_predecessors.md#checking-the-three-pairs), since $R-Q>Q-P$. It remains to exclude avoiding $r$ together with each of these primes. Fix $W\in\{P,Q,R\}$ and suppose a row of weight $dF+\epsilon$, $\epsilon=0$ or $1$, avoids both. In base $r$ the weight is $dr+dg+\epsilon$, with a units digit smaller than $P$. Write the parts $e_i r+j_i$. Carry freedom gives $\sum e_i=d$ and $\sum j_i=dg+\epsilon$. In base $W$ their low residues are nonnegative $\delta_i$ summing to $\epsilon$. As $r\equiv-g\bmod W$, we have
$$j_i\equiv ge_i+\delta_i\pmod W.$$
Both sides lie between zero and $Dg+1<W$, so they agree as integers. Thus all three parts are **exactly** $e_i F+\delta_i$. At most one positive part can have $e_i=0$ because $\sum\delta_i\le1$. At least two multipliers $e_i$ are positive. Dividing the parts above their units block $W$ shows that the quotients $e_iF/W$ must add without carry in the base prime of $W$.

* At $W=P$, for $e_i>0$ we have $e_i<B$ because there are at least two positive multipliers and $d\le B$. The low residues of $e_i QR$ modulo $B$ are $B-e_i$. Two such residues sum to at least $B$.
* At $W=Q$ use $PR=(Q-a)(Q+b-a)$. The low residues of $e_i PR$ modulo $Q$ are $Q-e_i a(b-a)$ for positive $e_i$. The two residues sum to more than $Q$ since $d b^2<Q$.
* At $W=R$ use $PQ=R^2-(2b-a)R+b(b-a)$. For $e_i>0$, the units coefficient $e_i b(b-a)$ is less than $R$, and the digit in the next block is $R-e_i(2b-a)$. Two such digits sum to more than $R$ because $d(2b-a)<R$.

Each case gives a carry. Consequently each trinomial row avoids at most one of the four fixed primes, so three rows share at least one of them.

## Minimum size for many multipliers in the example

Take $D=B=2^{75}$ in the displayed example. No prime divisor of $N=dF$ can be universal, where a universal prime divides every positive trinomial coefficient. In base 2 a weight has at most two digit units only if it is a power of two or of the form $2^i(1+2^j)$. The odd factor $Q$ of $N$ divides none of $2^j+1$ for $1\le j\le1+\lfloor\log_2(DF)\rfloor$; the verifier checks this finite assertion by modular doubling. Thus the binary digit sum exceeds two. For an odd prime $q\mid N$, a digit sum of at most two would require $N=q^i$, $2q^i$, or $q^i(1+q^j)$ with $j\ge1$. The large power $2^{256}\mid N$ leaves only the last form, with $j$ odd and $2^{256}\mid q+1$: an even power of an odd number is 1 modulo 4, and for odd $j$ the quotient $(1+q^j)/(1+q)$ is odd. The possible $q$ divide $d$, or equal $Q$ or $R$. All divisors of $d$ are less than $P-1$, and neither $Q$ nor $R$ is $-1\bmod P$. This excludes the last form too.

A prime $q\nmid N$ with digit sum at most two would require $N-1$ to be a power of $q$. The binary digit sum of the odd successor $N+1$ exceeds two, while in an odd base its digit sum is odd. Therefore a universal prime at the successor would require $N+1$ to be a prime power. The universality criterion follows by splitting base-prime digit units among three positive parts.

Require $3\parallel N-1$ and $5\parallel N+1$, or interchange 3 and 5. Here $\parallel$ means that the prime divides exactly once. Since both neighbors exceed five, these conditions exclude prime powers on both sides. As $F$ is a unit modulo 225, the two disjoint alternatives give $16$ multiplier classes modulo $225$ altogether. For each such $d$, every prime can be avoided by some row at each weight. Any list of three primes can then be defeated by selecting one row per prime. The four-prime list is minimum at both weights. There are exactly $2{,}686{,}501{,}821{,}365{,}842{,}610{,}459$ multipliers in the sixteen classes in the stated range. At $d=1$ the verifier also checks minimum size by independently enumerating all bases in which either weight could have digit sum at most two.

In this example $P\ge180Db^2$ also holds. The three pairs among the first three bases can therefore be checked by a second argument, taking $m=3$ in the [signed-polynomial test](signed_power_pairs.md#any-fixed-number-of-nearby-prime-powers).

Run `python3 check_fixed_four_primes.py` with python-flint installed. It proves primality of $Q,R,r$, checks the interval inequalities and congruences exactly, and checks the forced form of the parts against 23,377 exact smaller trinomial coefficients. The large interval is certified by the uniform proof. [FLINT primality test](https://flintlib.org/doc/fmpz/fmpz_is_prime.html)
