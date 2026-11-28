# Prime bases that work along whole geometric sequences of exponents

Fix a prime $p$. Write
$$
 \lambda_p=(2\sin(\pi/(4p+2)))^{-1},\qquad \beta_p=\log_p\lambda_p<1 .
$$
Let $I$ be any interval of $H\ge2$ consecutive positive integers. Let $\mathcal E$ be any nonempty finite set of positive integers coprime to $p$, and $Y\ge2$. Apart from at most
$$
 O_p\big(Y^2|\mathcal E|\,H^{\beta_p}\log(2HY)\big) \tag{1}
$$
prime bases $q\in I$, Wasserman's three-row conjecture holds **simultaneously** at all weights
$$
            1+D p^b q^{(p-1)e p^k},
 \qquad b\ge1,\ k\ge0,\ e\in\mathcal E,\quad
               D\le Y\text{ a prime power},\quad p\nmid D.
$$
The bound is uniform in the location of the interval and the sizes of the individual $e$. For one specified $D$, replace $Y^2$ in (1) by $D$, and $Y$ inside the logarithm by $D$.

For example, among primes $q\le X$ a proportion tending to one work at every
$$
      1+D\,2^b q^{e2^k},\qquad b\ge1,\ k\ge0,\quad
      e\le X^{1/5}\text{ odd},\quad D\le X^{1/20}\text{ an odd prime power}.
$$
Here $\beta_2=\log_2((1+\sqrt5)/2)<0.695$, so the count of possible exceptions is $O(X^{0.995}\log X)=o(X/\log X)$. Both unbounded exponents are allowed for the **same** base and every indicated cofactor. One may also replace “prime power” in both assertions by any of the cofactor types in the [four-prime extension](three_factor_exponent_rectangles.md#extension-to-four-predecessor-primes), taking $p,q$ coprime to $D$ for those additional cases.

## Two facts about the variable digits

First, for a prime base $q$ distinct from $p$ and coprime to $D$, a failure at $1+Dp^b q^j$ forces some $1\le x<D$ such that
$$
             xq^j+(D-x)q^j
$$
has no carry in base $p$. This is the [same-cut reduction](three_factor_exponent_rectangles.md#the-same-cut-in-both-bases). For the prime-power theorem coinciding bases are already covered by the result for predecessors with at most two primes.

Second, put $w(x)=(p/\lambda_p)^{v_p(x(D-x))}$. With **any fixed lower prefix** of a scale $t$, there are at most $w(x)\lambda_p^n$ choices for its next $n$ base-$p$ digits that can occur if $xt+(D-x)t$ adds freely in these positions. This is the [digit-block bound](three_factor_exponent_rectangles.md#how-many-exponent-terms-can-pass-a-cut); the proof counts just the indicated positions. When both slopes are units it counts at most $((p+1)/2)^n$ words. In particular the bound applies to the first $n$ positions even when higher positions carry. Also $\sum_{1\le x<D}w(x)\ll_p D$, by summing the valuation levels (at most $2D/p^a$ slopes have valuation at least $a\ge1$).

We use the elementary identity
$$
 v_p(u^j-v^j)=v_p(u-v)+v_p(j) \tag{2}
$$
for distinct positive units $u,v$ congruent modulo $p$ for odd $p$, or modulo 4 for $p=2$. In fact, raising to an exponent prime to $p$ leaves the valuation of the difference unchanged, by the geometric sum. Raising to $p$ increases it by exactly one, by binomial expansion and the extra factor of two in the binary hypothesis.

Fix any $j\ge1$ and $n=\lceil\log_p H\rceil$. Partition the prime-to-$p$ bases in $I$ into residue classes modulo $p^{f}$ with $f=1$ for odd $p$ and $f=2$ for $p=2$. In each class, the bottom $f+v_p(j)$ digits of $q^j$ are fixed. The next $n$ digits distinguish all the bases in that class by (2): a collision would require $p^{f+n}\mid q-q'$, impossible in an interval of length $H$. Applying the block bound and summing over the bounded number of classes and over $x$ gives $O_p(D H^{\beta_p})$ possible exceptional bases, simultaneously for all outer $b$. This estimate is independent of $j$.

## Replacing infinitely many exponents by one digit test

Fix $e\in\mathcal E$ and set
$$
        T_k(q)=q^{(p-1)e p^k}=1+p^{k+1}L_k(q).
$$
Fermat's little theorem and successive $p$th powers show $L_k(q)$ is an integer for every unit $q$. Expanding $T_{k+1}=T_k^p$ gives
$$
 L_{k+1}-L_k=\sum_{a=2}^p \binom pa
      p^{(a-1)(k+1)-1} L_k^a
      \equiv0\bmod p^k .
$$
Thus $L_k\bmod p^n$ is constant for all $k\ge n$, for each base. Choose $K=\max(n,\lceil\log_p D\rceil)$. For $k\ge K$, the low pieces $x,D-x$ in $xT_k,(D-x)T_k$ have sum $D<p^{k+1}$. Above the first $k+1$ places the cut is precisely $xL_k+(D-x)L_k$. Carry freedom at **any** such $k$ forces the first $n$ places of this addition to be carry-free at $L_K$. This implication uses only the fixed residue of $L_k$ and is uniform even if the cut $x$ changes with $k$.

Within each of the classes modulo $p^f$ above, the residues $L_K(q)\bmod p^n$ distinguish all bases in $I$. Indeed (2) gives $v_p(L_K(q)-L_K(q'))=v_p(q-q')-1$, since $p\nmid(p-1)e$. Equal residues would force $p^{n+1}\mid q-q'$. Applying the digit-block bound to the empty lower prefix therefore costs $O_p(D H^{\beta_p})$ bases for all indices in the tail $k\ge K$. Apply the fixed-exponent estimate at the $K$ earlier indices. Finally sum over $e$ and over eligible $D\le Y$, using $\sum_{D\le Y}D=O(Y^2)$. This proves (1).


## An explicit prime progression for all outer exponents

Let $M\ge2$, set $\epsilon=v_p(2)$ and $B=p^{M-\epsilon}$. Suppose the prime base satisfies
$$
                       q\equiv p^M-1\bmod p^M B .
$$
For any prime power $D$ coprime to $p$ and every $b\ge1$, the full conjecture holds at $1+Dp^b q^j$ under either condition:

- $j$ is odd and $D\le p^M$.
- $j$ is even and $D\cdot j/p^{v_p(j)}\le B$.

The same assertions hold for the additional cofactor types above when $p,q$ are coprime to $D$. For example, the prime $q=2^{61}-1$ works simultaneously at **every**
$$
                1+D\,2^b(2^{61}-1)^{2^k},\qquad
                b\ge1,\quad k\ge0,\quad D\le2^{60}\text{ an odd prime power}.
$$
In fact every odd exponent $j$ works at this base for $D\le2^{61}$, and every even exponent whose odd part is at most $2^{60}/D$ also works.

For odd $j$, set $t=q^j\equiv-1\bmod p^M$. The residues of $xt$ and $(D-x)t$ for any $0<x<D$ add to $2p^M-D\ge p^M$. This rules out the required cut.

For even $j$, write $q=-1+p^M u$, with $u\equiv1\bmod B$, and $j=p^v e$, $p\nmid e$. Since $j$ is even, $t=(1-p^M u)^j$. Its linear term gives
$$
        t=1+p^{M+v}L,\qquad L\equiv-e\bmod B. \tag{3}
$$
Indeed, $a\binom ja=j\binom{j-1}{a-1}$ gives $v_p(\binom ja)\ge v-v_p(a)$ for $2\le a\le j$. Every term of degree $a\ge2$ is divisible by $p^{2M+v-\epsilon}$, since $(a-2)M+\epsilon\ge v_p(a)$. The last inequality is equality for $a=2$ and follows from $v_p(a)\le a-2$ for $a\ge3$.

Both $x$ and $D-x$ are below $p^{M+v}$, so above these positions the parts of the cut are exactly $xL$ and $(D-x)L$. Modulo $B$ their residues sum to $2B-eD\ge B$, since both $ex$ and $e(D-x)$ are positive and strictly less than $B$. Thus the cut carries for every $x$. This proof works for arbitrarily large $v$.

Run **python3 check_prime_bases_at_geometric_exponents.py** for exact finite checks of stabilization, injectivity, and the common tail test, including carry tests from binomial coefficients. The verifier also confirms primality of $2^{61}-1$ by exhibiting 37 as a unit of order $2^{61}-2$ using the full factorization of that order. The proof of (1) uses the inequalities for all parameter sizes.
