# A restriction on the rows of a counterexample

For $N\ge3$ let $g=N-r$, where $r$ is the largest prime at most $N$. The following are necessary conditions on an individual selection of three rows with gcd one.

1. **A small part at every sufficiently large weight.** For all sufficiently large $N$, every counterexample has a part equal to 1 or 2 somewhere among its nine parts.
2. **A prescribed row at almost every weight.** If $N\ge12$ and a counterexample does not include $(1,1,N-2)$ up to permutation, then
   $$g\ \ge\ 120^{-7/5} N^{1/5}.$$
   In particular, for all but $o(X/(\log X)^A)$ integers $N\in[X,2X]$, for each fixed $A>0$, any counterexample must contain $(1,1,N-2)$. At those weights the entire conjecture reduces to the question whether **any two** positive trinomials have a common prime dividing $N(N-1)$.

The proof extends Bergman’s use of divisors of consecutive integers by counting the multiplicities of shifted factors of the parts. All small- and large-weight claims below are unconditional.

## A weighted divisibility lemma

Use $W(t)=t(t+1)/2$, and for integers $x\ge t\ge0$ define
$$B_t(x)=\prod_{j=0}^{t-1}(x-j)^{t-j}.$$
Thus $B_0=1$ and $1\le B_t(x)\le x^{W(t)}$ for positive $x$. Fix $1\le H\le N$ and put $L=W(H)$, $M=\operatorname{lcm}(1,\ldots,H)$. For each of the three rows $(x_1,x_2,x_3)$ choose its own three integer cutoffs $t_i\ge0$, with $\sum_i t_i=H$ and $t_i\le x_i$.

**Lemma.** If for each prime divisor $p$ of $\prod_{h=0}^{H-1}(N-h)$ some row is carry-free in base $p$, then
$$\prod_{h=0}^{H-1}(N-h)^{H-h}\quad\mid\quad M^L\prod_{\text{rows}}\prod_{i=1}^3 B_{t_i}(x_i).$$

**Proof.** Count the valuation one prime-power layer $p^e$ at a time. A layer contributes $\sum_{0\le h<H,\ p^e\mid N-h}(H-h)$ on the left. If $p^e\le H$, that contribution is at most $L$, covered by $M^L$. If $p^e>H$, at most one such $h$ occurs. In a row carry-free in base $p$, the residues $u_i=x_i\bmod p^e$ sum to $N\bmod p^e=h$. This follows by adding the lowest $e$ digits. If $u_i<t_i$, the factor $x_i-u_i$ contributes at least $t_i-u_i$ copies at this layer. The combined contribution of that row is at least
$$\sum_i\max(0,t_i-u_i)\ \ge\ H-h.$$
Summing over layers proves the divisibility. In particular, for a counterexample and $N\ge2H$,
$$N^L\ \le\ (2M)^L\prod_{\text{rows}}\prod_i B_{t_i}(x_i). \tag{1}$$

For a counterexample the lemma's hypothesis follows from Kummer. By Bertrand's postulate $r>N/2$ for $N\ge3$. A row which is carry-free in base $r$ has one part at least $r$. Designate one such row $(a,b,c)$, sorted increasingly. We will use $a+b\le g$. The other two rows need not be sorted unless a cutoff is assigned by size.

## Proof of the small-part restriction

Assume all nine parts are at least 3 and $N\ge36$. In estimates below the implied constants are absolute and come from (1) with $H=6$ or $H=18$.

If all parts are at least 6, take $H=18$ and all cutoffs equal to 6. Here $L=171$ and $W(6)=21$. Seven parts are at most $N$, and $a,b$ are at most $g$. Consequently $N^{171}\ll N^{147}g^{42}$ and $g\gg N^{4/7}$.

Otherwise let $m$ be the smallest part, so $m=3,4$ or 5. Take $H=6$, hence $L=21$. Give the two rows other than the row containing the selected occurrence of $m$ the cutoffs $(2,2,2)$. Each costs at most $N^9$ in (1). Give $m$ cutoff $m$, which supplies the absolute constant $B_m(m)$.

- For $m=4$ give the other two parts of that row cutoffs 1 and 1. For $m=5$ give them 1 and 0. The right side of (1) is $O(N^{20})$ in both cases, so $N$ is bounded.
- For $m=3$ give the other parts cutoffs 2 and 1. If the selected row is the designated row, give the other small part cutoff 2. Then $N^{21}\ll N^{19}g^3$. If the selected row differs from the designated row, its two remaining parts cost at most $N^4$, the designated row costs at most $N^3g^6$, and the last row at most $N^9$. Thus $N^{21}\ll N^{16}g^6$.

These alternatives give $g\gg N^{4/7}$ for every sufficiently large counterexample with all parts at least 3. On the other hand $g=O(N^{57/100})$. For example, apply Guth–Maynard’s Corollary 1.3 to an interval of length $N^{57/100}$ starting at $N-2N^{57/100}$, using $17/30<57/100$. Since $57/100<4/7$, this is a contradiction for sufficiently large $N$. [Guth–Maynard, Corollary 1.3](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

## Proof of the prescribed-row restriction

Suppose the three rows constitute a counterexample and none is $(1,1,N-2)$. Bergman’s Proposition 10 excludes $3\le a+b\le10$ under these hypotheses, and $a+b=2$ would be the excluded row. Thus $a+b\ge11$, in particular $b\ge6$. [Bergman, Proposition 10](https://arxiv.org/abs/0806.0607)

At most one row can contain a 1. In fact, for every prime $p\mid N$, a row avoiding $p$ must have all parts divisible by $p^{v_p(N)}$. If two rows contained a 1, the third would have all parts divisible by $N$.

Use $H=6$. Each of the two rows other than the designated row has a choice of cutoffs costing at most $N^9$. If it contains no 1, take $(2,2,2)$. If it contains a 1, order its parts increasingly and take $(1,2,3)$. Its second part is at least 2 by the excluded-row assumption, its third is at least 3, and $B_1(1)=1$ leaves total degree $3+6=9$ on the other two parts.

If $a=1$ take $(1,5,0)$ on the designated row. Its contribution is at most $g^{15}$. If $a\ge2$ take $(2,4,0)$, costing at most $g^{13}\le g^{15}$. With $M=60$ and $L=21$, (1) gives $N^{21}\le120^{21}N^{18}g^{15}$, as asserted.

For the frequency claim, let $J=\lfloor X^{0.18}\rfloor$. For large $X$, Guth–Maynard’s Corollary 1.4 supplies a prime in $(N-J,N]$ for all but $O(J+X\exp(-c(\log X)^{1/4}))$ weights in $[X,2X]$. The $O(J)$ term accommodates the lower boundary when shifting the interval start. Here $2/15<0.18<1/5$, so such weights satisfy $g<120^{-7/5}N^{1/5}$ for sufficiently large $X$. Finally, the coefficient of $(1,1,N-2)$ is $N(N-1)$, giving the stated two-row reduction. [Guth–Maynard, Corollary 1.4](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

## Reproducing the arithmetic check

Run `python3 check_weighted_lemma.py`. For weights 12 through 80 the script tests the large prime-power layers directly against exact binomial products, for permutations of the cutoff choices in the proof. It checks 54,452 instances, including 34,394 in which the required divisor exceeds the prime itself.
