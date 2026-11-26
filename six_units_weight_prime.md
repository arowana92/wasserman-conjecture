# Six digit units at a prime dividing the weight

Let $K_A=N!/\prod_i A_i!$ for a triple of positive parts of sum $N$. Write $s_b(N)$ for the base-$b$ digit sum, counting repeated digit units separately.

**A single coefficient.** If $b$ is prime, $b\mid N$, $s_b(N)\le6$ and $b\nmid K_A$, then
$$
\gcd\bigl(N(N-1),K_A\bigr)>2N.
$$
Consequently the full conjecture holds whenever two distinct primes dividing $N$ have digit sums at most five and six, respectively. If both digit sums are six, any counterexample must have a single row avoiding both primes.

**An even six-bit classification.** Suppose $N\equiv2\bmod4$ has exactly six ones in binary. If it admits a counterexample, it has the form
$$
N=(1+2^D)(2+2^j+2^k),\qquad 1<j<k<D/2 .
$$
Its row with odd coefficient consists of the three lower powers multiplied by $1+2^D$. In particular at most $O(L^3)$ weights in this class below $2^L$ can be exceptional, out of $\binom{L-2}{5}$ weights. The full conjecture holds for **every weight in the class below $2^{160}$**.

## Retaining divisors of the endpoint

Put $M=N-1$, $U=N/\gcd(N,K_A)$ and $S=M/\gcd(M,K_A)$. Then $U$ divides each part, since $K_A A_i/N$ is integral. At every full power in $S$ the residues of the three parts are $1,0,0$ in some order. To see this, take $r^e\parallel M$ and $v=v_r(K_A)<e$. If the sum of the three residues at $r^{e-v}$ exceeds 1, it continues to exceed 1 at every higher power through $r^e$. Each of these $v+1$ layers contributes at least one to Legendre's valuation formula, a contradiction. [Retained powers](endpoint_pair_bound.md)

Group the prime powers in $S$ by the position of the residue 1. Call their products the bins and omit bins equal to 1. We prove $2US<M$. If there is at most one bin, two parts are multiples of $US$ and the third is at least $U$. Since $b\nmid K_A$, $U$ contains the full power $B=b^{v_b(N)}\ge2$. Thus $N\ge2US+U$ gives the claim. In particular this argument applies whenever a part is a pure power of $b$, since it is coprime to $S$.

A row avoiding $b$ distributes its digit units among three nonempty parts without carries. Unless there is a pure power, the six units give exactly two to each part. Write these parts as $b^{l_i}F_i$ with $F_i=1+b^{d_i}$, $d_i\ge0$. Here $M$ is odd. Indeed an odd base contributes six odd units to $N$, and for $b=2$ we assumed $2\mid N$. If some $d_i=0$, its $F_i=2$ is a unit modulo $S$, which again leaves at most one bin.

Suppose there are at least two bins. A bin for which positions $i,j$ are zero divides $\gcd(F_i,F_j)$ and is coprime to $F_k$ for the remaining index $k$. All gap exponents are positive. Any odd prime dividing both $F_i,F_j$ forces $v_2(d_i)=v_2(d_j)$ by the order of $b$ at that prime. Two bins connect all three indices, so all have the same 2-adic valuation. Put $g=\gcd_i d_i$, $D=\max_i d_i$ and $g_{ij}=\gcd(d_i,d_j)$. Then
$$
\gcd(F_i,F_j)=1+b^{g_{ij}},\qquad \gcd_i F_i=1+b^g .
$$
For example, reduce modulo the numerical gcd and apply Bézout to the two exponents. Their quotients by $g_{ij}$ are odd, so an integer combination giving $g_{ij}$ has odd coefficient sum. This gives $b^{g_{ij}}\equiv-1$ and the reverse divisibility follows by factoring odd powers plus one.

The integer $1+b^g$ divides $N$ and is coprime to $S$. An occupied pair therefore has $t_{ij}=g_{ij}/g>1$, an odd integer, and its bin divides $(1+b^{g t_{ij}})/(1+b^g)<b^{g(t_{ij}-1)}$. The $t_{ij}$ for different pairs are coprime, since the gcd of the two corresponding pair gcds is $g$. The two largest occupied values $v<w$ belong to pairs sharing a part, so $D\ge gvw$. Also $\sum t_{ij}\le vw$: with two terms use $v,w\ge3,5$, and with three use $\sum t_{ij}\le2v+w<vw$ since then $v\ge5,w\ge7$.

Write $m$ for the number of bins. Since $U/B\mid1+b^g$,
$$
US < B(1+b^g)b^{g(\sum t_{ij}-m)}
\le \tfrac32 B b^{D-(m-1)g}.
$$
For $m=2$ the difference $vw-v-w\ge7$ gives at least seven additional powers of $b$ in the denominator of this last bound. For $m=3$, $(m-1)g\ge2$. In either case $2US<B b^D\le M$, as the largest-gap part is at least $B(1+b^D)$. This proves the coefficient theorem since $\gcd(NM,K_A)=NM/(US)$.

If distinct rows $A,C$ avoid two such sparse weight primes, their respective gcds with $NM$ exceed $2N$. Their common gcd with $NM$ exceeds $(2N)^2/(NM)>4$, since the lcm divides $NM$. By the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) a counterexample includes the endpoint of coefficient $NM$, so the two primes must both be avoided by one row. If either base has at most five units, its row contains a pure power of that base. Avoidance of the other weight prime would require all parts divisible by that other prime. This proves the assertions for two bases.

## Six binary units when the full power of two is 2

In a counterexample assign 2 to an avoiding row $A$ and assign all prime divisors of $NM$ to one of the two nonendpoint rows avoiding them. In this section $U,S$ denote the full factors assigned to $A$ from $N,M$. The complementary row has parts $(N/U)y_i$ with $\sum y_i=U$, so $U\ge3$. Each power in $T=M/S$ divides two of the $y_i$, giving
$$
M\le (U/3)^{3/2}S. \tag{1}
$$
Since $U=2u$ with $u>1$ odd, $A$ cannot contain a pure power of two. Its three parts are therefore $2^{l_i}(1+2^{d_i})$, and $u$ divides their three odd factors. This again gives equal 2-adic valuations for the $d_i$, now even if there are no bins in $S$. Thus $u\mid1+2^g$ and $(U/3)^{3/2}\le 2^{3g/2}$.

Two or more bins in $S$ are impossible. The bin bounds above give $M<2^{3g/2+g(\sum t_{ij}-m)}\le2^D$, whereas $M\ge2^{D+1}$. If there are no bins, (1) gives $M\le2^{3g/2}$, so $D=g$, as $D/g$ is odd. All gaps equal $D$ and $(2^{l_1}+2^{l_2}+2^{l_3})(1+2^D)-1\le2^{3D/2}$ yields $\max l_i<D/2$. The smallest $l_i$ is 1.

## Excluding exactly one bin

If only positions 1 and 2 are zero modulo $S>1$, write $d'=\gcd(d_1,d_2)$. The bin bound gives $S<2^{d'-g}$, so $2^{D+1}\le M<2^{d'+g/2}$ by (1). Since $D$ and $d'$ are multiples of $g$, $D=d'=d_1=d_2$. Now $D=gh$ and $d_3=gt$, where $h,t$ are coprime odd integers and $h\ge3$. The first two lower positions are distinct positive exponents, giving $M>6\cdot2^D$. In particular $2^{g/2}>6$ by the same upper bound, so $g\ge6$.

We have $S\mid 2^D+1$ and $S\mid 2^{l_3}(1+2^{gt})-1$. Choose a nearest integer $a$ to $l_3/g$, put $r=l_3-ag$ and $x=2^g$. The two polynomials
$$
F(X)=X^h+1,\qquad H(X)=2^{\max(r,0)}X^a(1+X^t)-2^{\max(-r,0)}
$$
both vanish at $x$ modulo $S$. Their integer resultant is a multiple of $S$ by the Bézout identity for the resultant. It is nonzero. At a root $\zeta$ of $F$, $(\zeta^t)^h=-1$. Hence $|1+\zeta^t|<2$ and $|1+\zeta^t|\ne1$, the latter because equality requires a primitive cube root. For $r<0$ or $r=0$ these facts exclude a zero of $H$. For $r>0$ a zero would make the algebraic integer $\zeta^a(1+\zeta^t)$ equal to the nonintegral rational $2^{-r}$.

Consequently $S\le(2\cdot2^{g/2}+1)^h$ by the product formula for the monic resultant. Write $y=2^{g/2}\ge8$. The ratio $(2y+1)/y^2$ is less than one, so the case $h=3$ suffices to show
$$
(2(1+y^2)/3)^{3/2}(2y+1)^h<6 y^{2h}.
$$
For $h=3$ divide by $y^6$. The left side becomes $(2/3)^{3/2}(1+y^{-2})^{3/2}(2+y^{-1})^3$, decreasing in $y$. At 8 its square is $(2/3)^3(65/64)^3(17/8)^6<36$. This contradicts (1) and the lower bound $M>6\cdot2^D$. The case without bins is the only one left.

## Checking the repeated blocks

For each possible gap and lower triple, let $F=1+2^D$. Assign to the row $(2F,2^j F,2^k F)$ every weight prime that it avoids. Moving such a prime to this row only relaxes the requirements on the companion. All predecessor prime factors must be avoided by the companion since each first-row entry is coprime to $M$. Test (1), then $2w+1\le U$ for every full power $w$ in $M$, using the two positive normalized entries divisible by $w$. Since $M>U$, each normalized entry $y$ is determined by its choices $y\equiv0,U\bmod w$. Check whether any three positive entries sum to $U$.

| Upper bound on $N$ (exclusive) | Repeated blocks passing the initial size bound | After actual assigned bases and (1) | Companion CRT systems | Admissible single entries | Triples summing to $U$ |
|---|---:|---:|---:|---:|---:|
| $2^{72}$ | 4,787 | 449 | 252 | 8 | 0 |
| $2^{96}$ | 12,527 | 1,009 | 620 | 12 | 0 |
| $2^{128}$ | 31,862 | 3,145 | 2,085 | 24 | 0 |
| $2^{160}$ | 64,874 | 6,354 | 4,162 | 30 | 0 |

A separate enumeration builds the actual part as $(M/w)((M/w)^{-1}\bmod w)$ for every unitary divisor $w$ of $M$, then checks disjoint groups and the scale of the actual row. This independent enumeration agrees on all four counts through $2^{160}$, covering 769,754,986 weights. It determines the scale by digit tests at all prime factors of the weight. The main checker uses factorial floors for the first row and CRT for normalized single entries. A third check directly tries all 395,010 pairings at the 26,334 six-bit weights below $2^{24}$ with $2\parallel N$. Of these, 59 survive the raw size test using gcds in the actual entries, and every survivor has equal gaps. The single-coefficient inequality was also compared with exact binomial products in 894 cases through weight 160.

~~~sh
# Python with python-flint installed
python3 check_six_even_pairs.py 160
python3 check_six_even_cross.py
~~~
