# Two translated factors with growing cofactors

Let $P=p^a$, $p$ prime, $a\ge1$, and let $0<s<h$ and $d\ge1$ be integers. Set $M=P(P+s)(P+h)$ and $L=2h+1$.

* **At the successor $dM+1$**, require just one prime power $w_h>1$ dividing $P+h$. Write $c_h=(P+h)/w_h$. The full conjecture holds if
$$P\ge\max\{c_h d^5 L^6,\ 16d c_h^3 L,\ 16d c_h^2 L^2\}. \tag{1}$$
There is no factorization condition on $P+s$.
* **At both $dM$ and $dM+1$**, choose prime powers $w_j>1$ dividing $P+j$ for $j=s,h$. Put $c_j=(P+j)/w_j$ and $C=\max(c_s,c_h)$. It suffices that
$$P\ge\max\{d^{16/3} c_s c_h L^6,\ 16d C^3 L,\ 16d C^2 L^2\}. \tag{2}$$
The prime bases need not be distinct. In particular the conclusion at both weights holds for either ordering of the gaps if the three factors themselves are prime powers and $P\ge d^{16/3}(2h+1)^6$.

For fixed shifts and $d$, these bounds allow cofactors as large as a sufficiently small constant times $P^{1/3}$. The [earlier translated-product bound](translated_products.md) allows $P^{1/6}$ for three factors. When the factors are prime powers and $h-s>s$, the [shared-prime theorem](clustered_predecessors.md) gives the larger range $2dh^2\le P$.

## The row avoiding the base of $P$

Use the [endpoint theorem](all_weights_endpoint_reduction.md). In a counterexample let $A$ receive the prime $p$, and let $C'$ be the other nonendpoint row. Write $UV=N$ and $ST=N-1$ for the full powers assigned to $A,C'$, in that order. Thus $U,V\ge3$ are coprime, and all parts of $A,C'$ are divisible by $U,V$, respectively.

Put $\epsilon=N-dM$, equal to 0 or 1. At either threshold the coefficients of $dX(X+s)(X+h)+\epsilon$ are smaller than $P$. Grouping the carry-free base-$p$ digits into blocks shows that
$$A_i=P F_i(P)+\delta_i,\qquad F_i(X)=v_iX^2+a_iX+b_i,$$
where the nonnegative integers $v_i,a_i,b_i$ sum over $i$ to $d,d(s+h),dsh$, respectively. The nonnegative $\delta_i$ sum to $\epsilon$. Each $F_i$ is nonzero by positivity and $U\ge3$.

Both cubics $G_\pm(X)=dX(X+s)(X+h)\pm1$ are irreducible over the rationals. Indeed a nontrivial integer factorization into $J,K$ would give $J=\pm K$ with the same sign at each of $0,-s,-h$, because the product takes the same value $1$ or $-1$ at all three nodes. The difference or sum has degree at most two and would vanish identically, contrary to odd degree. The cubics are primitive, so Gauss's lemma applies. Each complex root $z$ has $|z|\le h+1$ because $d|z(z+s)(z+h)|=1$.

Consequently for every part polynomial,
$$1\le |\operatorname{Res}(G_\pm,F_i)|\le d^5 L^6. \tag{3}$$
Indeed $|F_i(z)|\le d(h+1)^2+2dh(h+1)+dh^2=d L^2$. Multiply over the three roots and include $d^{\deg F_i}\le d^2$. The gcd of the two polynomial values at $P$ divides the resultant by the integral Bézout identity.

## The successor needs only the largest translated factor

Use $N=dM+1$. For either unshifted part ($\delta_i=0$), $U$ divides $F_i(P)$ since $\gcd(U,P)=1$. Also $U\mid G_+(P)$, so (3) gives $U\le d^5L^6$. Under (1), $w_h>U$. Row $C'$ cannot avoid its base: at least two normalized parts of $C'/V$ would be divisible by $w_h$, whereas their total is $U$. Thus $A$ must avoid this base too.

Write $h-s=t$, $c=c_h$, $w=w_h$. Substitute $X=cZ-h$ in each $X F_i(X)+\delta_i$. In degrees $0,1,2$ the resulting signed coefficients have absolute values at most $4dh^3+1,\ 8dch^2,\ 5dc^2h$, respectively, by direct expansion. Each is less than $w/4$ under (1). The total polynomial in the new variable is
$$1+dcZ(cZ-h)(cZ-t).$$
Its constant digit in base $w$ is 1 and its next digit is $dcht<w/4$. Carry freedom at the base of $w$ therefore forces each constant coefficient to be 0 or 1. A negative constant would have a residue greater than $3w/4$. It then forces every linear coefficient to be nonnegative by the same argument at the next digit, now without an incoming borrow.

Each of the two unshifted polynomials is nonzero with positive leading coefficient and has a positive real root $Z=h/c$ since its original constant term is zero. It must have a negative coefficient. Only the quadratic coefficient could be negative. Thus two different parts have a digit greater than $3w/4$ at the quadratic position. This contradicts carry freedom. This proves (1).

## The unshifted weight

At $n=dM$, let $w=w_j$, $c=c_j$ for either $j=s$ or $j=h$. A positive **trinomial** row cannot avoid both the base of $P$ and the base of $w$ under the last two bounds of (2). We give the argument also when the bases coincide. Avoidance forces the parts to be $A_i=P w f_i$. For a shared base, the full power dividing $n$ contains the product $Pw$ since $w\mid P+j$.

Put $u_i=\lfloor A_i/P^3\rfloor$. Their sum is $d$ and the remainders sum to $H=d(M-P^3)\le3dhP^2<w^3$. Also $(P+j)^3-M\le(P+h)^3-P^3\le12hP^2$, so $\lfloor n/w^3\rfloor\ge dc^3-1$. Here $w\ge P/c$ and $P\ge16dc^3 L$. At every position we have $\lfloor A_i/w^3\rfloor\le c^3u_i$.

Define $e_i=f_i-c^2w u_i$. Since $P w(c^2w)=P(P+j)^2>P^3$, the sum of the positive $e_i$ is at most $H/(Pw)$. In particular $j\sum_i\max(0,e_i)\le 3dc h^2<w$, by the last bound of (2). After dividing the row by $w$, its low residues at $w$ are $[-j e_i]_w$, since $P=cw-j$. Two positive $e_i$ would already contribute at least $w$ to their sum. Hence at least two of the three $e_i$ are nonpositive. At each such position, $u_i\ge1$ and
$$A_i\le u_iP(P+j)^2<u_i(P+j)^3=c^3 w^3 u_i.$$
The sum of floors at $w^3$ is now at most $dc^3-2$, contradicting the required floor of $n$.

In a counterexample this puts both chosen powers on $C'$, so $w_s w_h\mid V$ and $U\le dP c_s c_h$. If they share a base, its power in the product $n$ includes both contributions. If a base coincides with $p$, the preceding incompatibility has already excluded the counterexample.

Each full power in $S$ divides at least two of the $F_i(P)$, since it divides two actual parts and is coprime to $P$. With $G_-(P)=n-1$, (3) yields
$$S^2\mid\prod_i\gcd(n-1,F_i(P)),\qquad S^2\le d^{15}L^{18}.$$
Likewise $T^2$ divides the product of the three normalized companion parts of sum $U$, so $27T^2\le U^3$. Combining and using $n-1>dP^3$ gives $27P^3<d^{16}(c_s c_h)^3L^{18}$. This contradicts (2), which also implies (1).

## An example beyond the sixth-root cofactor range

Take $s=1$, $h=2$, and
$$P=1000000000046111699137200127,\qquad
P+1=2^{20}\cdot953674316450225543153,\qquad
P+2=3^{12}\cdot1881676423245688042769.$$
The three displayed numbers other than the cofactors are prime. Both cofactors exceed $P^{1/6}$. Inequality (2) certifies both neighbors for every $1\le d\le118$. At $d=119$ the first inequality in (2) no longer holds. This is a limit of the stated sufficient condition.

Run `python3 check_translated_cubic_bounds.py` with python-flint installed. It proves the three primalities by recursive modular order certificates and checks the inequalities at all 118 multipliers. It also verifies 333,452 resultant bounds by exact determinants. At smaller scales it checks 1,117 incompatible assignments by enumerating candidate rows, finds 22 compatible controls outside the sufficient conditions, and verifies the predicted carry at the quadratic digit in 162 shifted rows.
