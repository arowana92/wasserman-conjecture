# Separation of binary exponents at cofactor five

Fix any odd prime power $Q$ and write $N_k=1+5\cdot 2^kQ$ for $k\ge1$. **Among any 1,356 consecutive binary exponents, the full three-row conjecture can fail at most once. Among any 4,001 consecutive exponents it can fail at most twice.** Both bounds are uniform over $Q$ and the location of the interval.

More precisely, if $s$ is the positive difference of two failing exponents for the same $Q$ and $s\le4000$, the only possibilities left here are
$$
                       1356,\quad 1424,\quad 2936,\quad 3364.
$$
No two of these numbers sum to a member of the list, even allowing repetition. Combining the first bound with the [separation of odd powers](five_power_separation.md) gives at most two possible failing pairs in the entire array of any 1,356 consecutive positive binary exponents and the first 1,000 odd powers of any given prime power $R$.

There is also a restriction on arbitrarily large differences. The set of positive $s$ that could occur for **any** $Q$ has upper natural density **less than 6 percent**. Inside every fixed arithmetic progression of positive differences, its relative upper density is **less than 61 percent**.

## A binomial theorem behind the separation

Put $B(t)=\binom{5t+1}{3t}$ for integers $t>0$. The same finite list and both density bounds apply to
$$
  \mathcal D=\{s\ge1:\text{ for some }t\ge1,\quad 5\nmid B(t)B(2^s t)\}.
$$
Indeed a failure at $N_k$ requires a trinomial coefficient coprime to 5 with one part $3t$, $t=2^kQ$, by the [cofactor-five reduction](predecessor_five_seven.md#full-factor-five). Factoring off the binomial $B(t)$ shows that $5\nmid B(t)$.

By Kummer, $5\nmid B(t)$ precisely when $3t$ and $2t+1$ add without a carry in base five. If $5\mid t$, their digits at the first nonzero place of $t$ sum to 5. If $t$ is a unit, checking the units digit forces $t\equiv2\pmod5$. Thus $s\in\mathcal D$ implies $4\mid s$. Kummer's criterion for the original problem also appears in [Bergman's paper](https://arxiv.org/abs/0806.0607).

## Low and leading digits

For $D=2^s$, first extend the base-five residue $t\equiv2\bmod5$ one digit at a time. At modulus $m=5^j$, keep the extension only if
$$
                [3u]_m+[2u+1]_m<m\qquad(u=t,Dt),
$$
where brackets give the least residue. If a level is empty, the exclusion holds at every scale. For instance $s=36$ first has an empty level at 262. The four lowest places exclude every $s\equiv4\bmod500$, since $2^{500}\equiv1\bmod625$. The finite implication on the 625 residues and its extension to every exponent in this progression were checked in [Lean](FiveExponentGaps.lean).

For the leading test, suppose $5\nmid B(t)$ and put $u=t+1/2$. For every $m=5^j$ with $j\ge1$,
$$
    \left\lfloor\frac{5u}{m}\right\rfloor
     =\left\lfloor\frac{3u}{m}\right\rfloor+
      \left\lfloor\frac{2u}{m}\right\rfloor .
$$
In fact $5t\bmod m$ is a multiple of 5, so increasing $5t$ by $1$ or $5/2$ preserves its quotient. On the right $2u=2t+1$, and increasing $3t$ by $3/2$ can only raise its quotient. The carry-free equality before this increase and the floor inequality afterwards prove the identity. For $j\le0$, the fractional part of $u/5^j$ is $1/2$, which also satisfies the identity. Checking the breakpoints of $\lfloor5y\rfloor-\lfloor3y\rfloor-\lfloor2y\rfloor$ now gives
$$
 \{u/5^j\}\in A=[0,1/5]\cup[1/3,2/5]\cup[1/2,3/5]\cup[2/3,4/5]
                   \qquad (j\in\mathbb Z).
$$
The closed endpoints enlarge the allowed set.

Normalize $x=(t+1/2)/5^n\in[1/5,1)$, and choose $b$ with $5^b\le D<5^{b+1}$. At each level $i$, the first mantissa $5^ix$ must belong to $\mathbb Z+A$. For the second binomial apply the same test to
$$
   5^i(D/5^b)(x-\delta),\qquad
   \delta=(1-1/D)/(2\cdot5^n).
$$
For a chosen height $h$ and $t\ge5^{h+7}$ grid the $x$ axis in units $1/(30D5^h)$. The endpoints of both unperturbed tests are grid integers for $-1\le i\le h$. In these coordinates $0\le\delta\cdot30D5^h\le15D/5^7$. Enlarge the right endpoint of every second interval by $\lceil15D/5^7\rceil$ grid units. An empty intersection certifies all such $t$. Below $5^{h+7}$ enumerate every allowed low prefix of $h+7$ digits and check the remaining higher digits of $3t,2t+1$ and their counterparts at $Dt$.

Of the 1,000 multiples of four through 4,000, the low test excludes 943 and the leading test with its finite lower range excludes 53 others. At $s=464$ it takes 296 leading levels above level zero. All 407 possible short inputs fail the complete binomial test. The other four differences are undetermined by this calculation.

## Density of possible differences

At $m=5^9$ there are 1,145 unit residues passing all nine low layers. Write them $2\cdot16^j\bmod m$. Since $v_5(16-1)=1$, successive fifth powers of 16 increase this valuation by one, and the order is $5^8$. The pairwise differences between the exponents $j$ occupy exactly 151,047 of the $5^8=390,625$ classes. Every $s\in\mathcal D$ must be four times one of those classes modulo $4\cdot5^8$.

For the real restriction, let $C$ be the union obtained by intersecting $[1/5,1]$ with the seven conditions $5^ix\in\mathbb Z+A$, $0\le i\le6$. For each pair of component intervals $[a,b],[c,d]$, take $[a/d,b/c]$ and its fivefold multiple, both clipped to $[1,5]$. This represents the ratio of two mantissas modulo multiplication by five. The endpoints 1 and 5 represent the same point on the logarithmic circle. The union omits more than $0.39\log5$ in logarithmic length. To certify this inequality, on each omitted interval $(v,w)$ use $\log(w/v)\ge(w-v)/w$. Even after rounding the included intervals outwards on a $2^{-36}$ grid and rounding each lower bound down to that grid, the sum is at least $43273786538/68719476736>0.39\cdot1.61$. The inequality $\log5<1.61$ follows from the first ten terms of the exponential series.

For a witnessing $t$, the logarithm of $(Dt+1/2)/(t+1/2)$ differs from $\log D$ by at most $\log(1+1/(2t))$. Hence witnesses above any cutoff have $s\log_5 2$ modulo one in a fixed neighborhood of the closed ratio set on the logarithmic circle, with radius tending to zero as the cutoff increases.

Each fixed smaller $t$ contributes a zero-density set of $s$, even in any fixed arithmetic progression. In a successful input $Dt$, a digit 4 in any place would put the fractional part of $(Dt+1/2)/5^j$ strictly between $4/5$ and 1 at that place. There are at most $4^j$ allowed residues modulo $5^j$. For a unit $t$, the powers $2^s t$ cycle through $4\cdot5^{j-1}$ residues, and using a fixed progression of step $v$ reduces the period by a factor at most $v$. Taking $j$ to infinity proves the zero-density claim. Nonunits cannot witness a gap.

Finally $\log_5 2$ is irrational, so the fractional parts of $s\log_5 2$ are equidistributed on every fixed progression. This follows by summing a geometric series for each nonzero Fourier frequency. First take the upper density at a fixed cutoff and then shrink the neighborhood. The finite union of closed arcs has relative length strictly below $0.61$. This proves the progression bound, including on each of the allowed congruence classes from the low test. Globally the upper density is at most $0.61\cdot151047/(4\cdot390625)<0.06$.

## Reproduction

The first program generates the proof heights for each difference. The independent verifier checks all 1,000 records. It reads base-five digits by Lucas' theorem, tracking the multiplication by $D$ and comparing the digit of $3t$ with the preceding input digit in $5t+1$. Its leading calculation uses $t/5^n$ directly. Within a cell $\lfloor5y\rfloor=k$ it bounds $\lfloor3y\rfloor+\lfloor2y+\epsilon\rfloor$ by solving the inequalities $3y\ge j$, $2y+\epsilon\ge k-j$ for integer $j$, using the maximum possible offset from the added 1. It independently computes the finite residue ratios by multiplication and the logarithmic gaps with exact rational endpoints.

~~~sh
mkdir -p .work
python3 check_five_exponent_gaps.py .work/five_gap_proofs.json
python3 check_five_exponent_gaps_cross.py .work/five_gap_proofs.json
lean FiveExponentGaps.lean
~~~
