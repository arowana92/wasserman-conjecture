# Successors of products with a nearby square or higher prime power

**Theorem.** Let $p$ be prime, $a\ge1$, $j\ge2$ and $1\le d\le1000$. At
$$N=1+p^a(p^j+d),$$
Wasserman's three-row conjecture holds for all of these parameters except possibly when $(p,j,d)=(37,2,909)$. Thus every $j\ge3$ works throughout the offset range, and every $j\ge2$ works through offset 908. Both exponents are unrestricted in these conclusions.

On the remaining family $N=1+2278\cdot37^a$, the conjecture holds for every
$$a\le178666130578561716682765137216.$$
Any exceptional exponent would satisfy $a\equiv2172940\pmod{2718912}$. Any two exceptional exponents differ by at least $94{,}438{,}145{,}595{,}691{,}968$.

## Exhausting the row assignments

Put $P=p^a$, $H=p^j$, and $Q=H+d$. By the [endpoint reduction](all_weights_endpoint_reduction.md), it suffices to treat the endpoint $(1,1,N-2)$ and two other rows $A,C$. If the three coefficients have gcd one, assign each prime of $N(N-1)$ to a row among $A,C$ avoiding it, with $p$ assigned to $A$. Write $UV=N$, $PBT=N-1$, where $U,PB$ are the products of the full powers assigned to $A$ in $N,N-1$, respectively. Then $B T=Q$ and $C=V(f_1,f_2,f_3)$ with $f_1+f_2+f_3=U$. All three parts of $A$ are multiples of $U$. In particular $U\ge3$ and $\gcd(U,Q P)=1$.

The [successor size bound](successor_powers.md) gives $27Q^2<d^5$. For completeness, the base-$p$ digits force $A=(\alpha P,\beta P,\gamma P+1)$ with positive slopes adding to $Q$ without a carry. Some slope receives a digit of place value at least $H$, even when $d\ge H$. Designate it large. The other two slopes $x,y$ satisfy $x+y\le d$. Choose $y$ unshifted. Of the full prime powers of $B$, let $e$ be the product of those whose exceptional residue 1 occurs at $Py$ and $f$ the product of the rest. The $p$-part of $Q$, if present, lies in $f$. Thus $e\mid Py-1$, $Uf\mid y$, and $e,f,T$ are pairwise coprime. If the large part is shifted, $Ue\mid x$. If $x$ is shifted, the large unshifted slope $Q-x-y$ is divisible by $Ue$, so $x+y\equiv0\pmod e$ and $x+y\equiv Q\pmod U$. Either way $e\le d$.

Every full power $w$ in $T$ divides two of the positive $f_i$. Consequently $2w+1\le U$ and $T^2\mid f_1f_2f_3\le U^3/27$. Combining these bounds gives $Q=efT\le dy\sqrt{U/27}<d^{5/2}/\sqrt{27}$. There is a sharper bound when the large part is shifted. Then $Ue\mid x$ and $Uf\mid y$, hence
$$Q\le\frac{xy}{\sqrt{27U}}\le\frac{d^2}{4\sqrt{27U}}\le\frac{d^2}{36}.$$

Enumerate $H=p^j$ starting at $p^2$ and $Q=H+d\le\lfloor\sqrt{1000^5/27}\rfloor$. Each full power of $Q$ is at most $d$ by the bounds on $e,f,T$. We omit the [settled cofactors of full prime powers](prime_power_cofactors.md) and weights whose predecessor has at most two distinct prime divisors. For every other $(p,H,d)$ assign each full power of $Q$ to $e,f$ or $T$, placing any $p$-part in $f$. Try all $U$ satisfying the size and coprimality bounds, all $y\in Uf\mathbb Z$ with $0<y<d$, and both progressions for $x$ above. When both small slopes are unshifted we may take $y\le x$.

For each resulting shape, impose $QP+1\equiv0\pmod U$ and the required residues $0,0,1$ at the powers assigned to $A$ in $B$. This gives $P\equiv R\pmod M$, with $U\mid M$. Intersect with the powers of $p$ modulo $M$ and test the possible low digits at each required prime of $A$. The normalized companion has positive parts of sum $U$ and, for each full $w$ in $T$, two parts divisible by $w$. Enumerate it by both exceptional positions and individual parts. The required bases of $C$ include the primes of $T$, and also 2 when $p,Q,U$ are all odd because then $2\mid V$. The [simultaneous exponent test](successor_exponents.md) uses the generalized Chinese remainder theorem on the possible exponent classes at all required bases.

| Stage | Entries |
|---|---:|
| First-row assignments lying in the power orbit | 32,799 |
| First-row assignments passing the separate low-digit tests with modulus cap 512 | 11,584 |
| Pairs after enumerating the normalized companion | 1,278 |
| Pairs passing simultaneous tests at modulus caps 256 and 4,096 | 32 |

At a cap, use the largest power of the base at most the cap, or the base itself if larger. At $r\mid U$ enforce $r\nmid V$ whenever the available precision determines it.

## Testing whether a digit prefix can finish

At a required base $r\ne p$ and $w=r^k$, enumerate $v=p^a\bmod w$ compatible with the class modulo $M$. The scalar for $A$ is $v$, and the scalar for $C$ is $(Qv+1)U^{-1}\bmod w$ whenever $r$ is assigned to $C$. Call the relevant scalar $z\in[0,w)$, write the required row as $s_i z+b_i$, and test its $k$ low digits. If those digits are carry-free, the multiplication carries entering the next column are $c_i=\lfloor(s_i z+b_i)/w\rfloor$. For an unknown next scalar digit $t$, the transition is allowed if
$$\sum_i ((s_i t+c_i)\bmod r)<r,$$
and the successor carries are $\lfloor(s_i t+c_i)/r\rfloor$. The carries stay bounded by the larger of the corresponding slope and initial carry. A finite high suffix exists exactly when some reachable state has carry-free entries, since the remaining scalar digits may then be zero. A real counterexample must pass this test at each required base.

At cap 4,096 this excludes 29 pairs, with 18 exclusions at 2, one at 3, one at 7, two at 41, six at 83 and one at 101. At $Q=5447,p=67$, the companion is $V(1,1,1)$ with $U=3$. Here $2\mid V$ and its first nonzero binary digit forces a carry. At $Q=38141,p=193$, the companion is $V(1,2,4)$ with $U=7$. Here $3\mid V$ and for each nonzero digit $t$ in base 3, the residues of $t,2t,4t$ sum to at least 3. These two contradictions eliminate all but one pair.

## The remaining exponent class and the finite interval

At $(p,H,d)=(37,1369,909)$ the only pair, up to reordering, has
$$A=(833P,1407P,38P+1),\qquad C=V(1,2,4),\qquad V=(2278P+1)/7.$$
Its progression $P\equiv5665\pmod{7973}$ gives $a\equiv28\pmod{48}$. The first row must avoid $7,17,67$ and the companion must avoid 2. Read the indicated low digits, retaining just prefixes with a finite carry-free completion at that base. Exhausting one exponent period gives the following table.

| Base | Digits | Period including the initial class | Admissible exponents modulo that period |
|---:|---:|---:|:---|
| 2 | 9 | 384 | 76, 268 |
| 7 | 3 | 2,352 | 2,044 |
| 17 | 3 | 13,872 | 8,908 |

Write $b=2172940$, $L=2718912$. The first condition is $a\equiv76\pmod{192}$, and combining the three gives $a\equiv b\pmod L$.

Lift this class through 30 base-seven digits, discarding incomplete prefixes at each step. Since $v_7(37^L-1)=3$, at precision $7^k$, $k\ge3$, each surviving class of period $S=L7^{k-3}$ has seven lifts. Replacing $a$ by $a+iS$ leaves its first $k$ scalar digits fixed and adds $i\,(37^a\bmod7)\,(37^S-1)/7^k$ to its next digit modulo 7. The increment for $i=1$ is a unit by binomial expansion. Each allowed next digit therefore determines exactly one lift.

This yields 878,080 distinct classes modulo $L7^{27}=178666130578561716682765137216$. Sorting them and including the gap across the period boundary gives the stated minimum separation. Take the least nonnegative representative of each class. Only 70,603 pass the first 18 binary digits of the companion. Every one of those representatives has a carry within the first 100 digits in at least one required base. The period itself is outside the initial exponent class. This proves the interval claim without constructing the full weight.

## Why a fixed number of digits cannot exclude the remaining family

Here is a simultaneous lifting criterion. Fix an integer $p>1$, a progression $a\equiv a_0\pmod{L_0}$ with $L_0\ge1$, and a finite nonempty set $S$ of primes coprime to $p$. For each $r\in S$ choose a set $W_r$ of $r$-adic integers as possible local witnesses for $p^a$. Let $L$ be any multiple of $L_0$ and of all the orders of $p$ modulo $r$ for odd $r\in S$, including its order modulo 4 if $2\in S$. Set $h_r=v_r(p^L-1)$.

For every prescribed precision $k$ there exist arbitrarily large positive $a$ in the progression with $p^a$ congruent modulo $r^k$ to some element of $W_r$ for every $r\in S$ **if and only if** there is a single $b\bmod L$ in the progression and, for each $r$, an element $w_r\in W_r$ with
$$w_r\equiv p^b\pmod{r^{h_r}}.$$
For necessity take $k\ge\max_r h_r$ and reduce one such exponent modulo $L$. For sufficiency fix witnesses at these initial precisions. For $t\ge0$, binomial expansion gives $v_r(p^{Lr^t}-1)=h_r+t$, using $h_2\ge2$ at 2. Thus varying $n\bmod r^{\max(0,k-h_r)}$ in $p^{b+Ln}$ gives every lift of $p^b$ modulo $r^k$ compatible with its initial digits. These conditions on $n$ have pairwise coprime moduli, so the Chinese remainder theorem supplies arbitrarily large simultaneous exponents. In particular a fixed choice of local witnesses works for all requested precisions.

Apply the criterion with $p=37$, the $b,L$ above, and the integer witnesses below. Each scalar satisfies the entire shape congruence $P\equiv5665\pmod{7973}$ and $49\nmid2278P+1$. At each odd base the listed scalar gives a carry-free finite first row. At 2 the companion scale is $V=144432397331726913$, and the binary expansions of $V,2V,4V$ are disjoint.

| $r$ | $r^{h_r}$ | $37^b\bmod r^{h_r}$ | Integer witness for $P$ |
|---:|---:|---:|---:|
| 2 | 256 | 209 | 443822116471505 |
| 7 | 343 | 226 | 1676335771086440963 |
| 17 | 4,913 | 905 | 27225360927 |
| 67 | 67 | 37 | 11694083 |

At 37 the first row always adds without carries by its normalized digits. Thus arbitrarily large exponents pass **any prescribed finite number of low digits** of the assigned rows simultaneously at the four other fixed bases. Each prefix can even be completed separately to an ordinary integer scalar satisfying the shape congruence and the full local carry test. The verifier constructs one simultaneous exponent of 417 decimal digits for a precision of 100 digits at every base.

### Excluding every auxiliary prime through ten million

For each prescribed digit depth the exponents above can also be chosen so that **every prime factor of $V=(2278\cdot37^a+1)/7$ exceeds $10^7$**. In particular no uniform finite digit cutoff at primes through ten million can rule out this shape by Kummer's criterion. The assertion concerns low digits. A common prime at any particular weight can still arise from a higher carry.

Write $a=b+Ln$ and $S=\{2,7,17,67\}$. For $r\in S$ take $m_r$ to be the least power of $r$ at least $10^7$. For the scalar $z_r$ in the table above, the following classes solve $37^{b+Ln}\equiv z_r\pmod{m_r r^{h_r}}$. The order calculation proves uniqueness modulo $m_r$.

| $r$ | $m_r$ | $n\bmod m_r$ |
|---:|---:|---:|
| 2 | 16,777,216 | 16,730,577 |
| 7 | 40,353,607 | 2,128,588 |
| 17 | 24,137,569 | 12,821,839 |
| 67 | 20,151,121 | 297,731 |

Let $m=\prod_{r\in S}m_r$ and $n_*$ be the simultaneous residue modulo $m$. At every prime $q\le10^7$ outside $S\cup\{37\}$ put $g=\gcd(q-1,m)$ and $f=(q-1)/g$. Then $f$ is coprime to $m$. Require $n\equiv0\pmod f$ at all these primes, which is compatible with arbitrary further lifts at the four primes of $S$. The two conditions $t\equiv n_*\pmod g$, $t\equiv0\pmod f$ determine $n\bmod(q-1)$. Direct modular exponentiation gives
$$2278\cdot 37^{b+Lt}+1\not\equiv0\pmod q$$
at all 664,574 auxiliary primes. Three of the four primes of $S$ and 37 divide $N-1$, and $7$ divides $N$ exactly once. This proves the assertion. The [independent checker](check_higher_successors_cross.py) tests primality directly, assembles the exponent from its separate CRT components, and verifies the first 5,000 by a single gcd with their primorial.

## Reproduction

The [Lean proof](HigherSuccessorBounds.lean) checks the deductions $27Q^2<d^5$ and $36Q\le d^2$ from their respective factor and part inequalities. Run it with Lean 4.36.0 using `lean HigherSuccessorBounds.lean`. The enumeration uses Python 3.9 or later and python-flint. Run the following from the repository root.

```sh
python3 check_higher_successors.py --list-local
python3 check_higher_successors_cross.py
python3 check_higher_successor_ray.py
```

The cross-check enumerates the small parts directly and tests the actual residues throughout the scalar orbit at each of the 23 assignments reaching the final stage. It compares the entire set of possible scale residues for each of the 65 first rows there. It also checks all 29 prefix exclusions using reverse reachability from accepting carry states. The exponent verifier compares the first five new lifting levels with direct modular exponentiation and checks each local witness by factorial floors.
