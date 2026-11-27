# Common exponent divisors at multiplier seven

**Let $p,q$ be any odd primes and $a,b$ any positive integers. If $a$ and $b$ have a common prime divisor at most 47, Wasserman's three-row conjecture holds at $1+7p^a q^b$.** For example this covers $1+7(pq)^{47}$ for every pair of odd primes, as well as $1+7p^{5u}q^{5v}$ for arbitrary positive $u,v$.

**For a set of primes $p$ of relative natural density one, the conclusion holds simultaneously for every odd prime $q$, every $a\ge1$, and every $b$ divisible by any prime at most 47.** These $b$ have natural density $9968041656757/11573306655157$, about $86.13\%$. The range can grow with $p$ when the two exponents have different prime factors, as detailed below.

## Perfect powers and sums of six repunits

Write $t=p^a q^b$. The [even-weight reduction](predecessor_five_seven.md#a-restriction-on-every-base-seven-digit) gives necessary conditions for failure:
$$t\equiv11\pmod{72},\qquad 6(t+1)=\sum_{i=1}^6 7^{e_i},\qquad e_i\ge1.$$
If a prime $\ell$ divides both exponents then $t=y^\ell$. Squares are excluded modulo 8. Cubes of integers prime to 3 are $\pm1$ modulo 9, whereas $t\equiv2\pmod9$. For each prime $5\le\ell\le47$ we prove the stronger arithmetic statement that there is no nonnegative integer $y$ with $y^\ell\equiv11\pmod{72}$ and $6(y^\ell+1)$ a sum of six powers of seven, even allowing zero exponents and repetitions.

At a period $h$, fold the exponents to $u_i=e_i\bmod h$ and put $T=\sum_i(7^{u_i}-1)/6$. Necessarily
$$T\equiv11\pmod{\gcd(432,7^h-1)/6}.$$
Indeed $t-T$ is divisible by $(7^h-1)/6$. For every prime $r\equiv1\pmod\ell$ dividing $7^h-1$ with $r\nmid6$, the residue $T\bmod r$ must be an $\ell$th power, including zero. For nonzero residues this is equivalent to $T^{(r-1)/\ell}\equiv1\pmod r$.

Start with every unordered multiset of six residues modulo $\ell$. At successive periods $2\ell,6\ell,12\ell,36\ell$, lift each surviving multiset by all choices of representatives of its elements. Apply the congruence and power tests with the primes listed in the verifier. Every solution must give a surviving multiset at each period. The exact counts follow.

| $\ell$ | $h=\ell$ | $2\ell$ | $6\ell$ | $12\ell$ | $36\ell$ |
|---:|---:|---:|---:|---:|---:|
| 5 | 45 | 22 | 38 | 5 | 0 |
| 7 | 19 | 3 | 0 | | |
| 11 | 76 | 13 | 1 | 0 | |
| 13 | 1456 | 170 | 0 | | |
| 17 | 257 | 181 | 0 | | |
| 19 | 421 | 17 | 0 | | |
| 23 | 45 | 37 | 17 | 0 | |
| 29 | 87 | 36 | 0 | | |
| 31 | 69 | 1 | 0 | | |
| 37 | 136 | 5 | 1 | 3 | 0 |
| 41 | 184 | 1 | 0 | | |
| 43 | 6832 | 8 | 40 | 0 | |
| 47 | 9114 | 3231 | 0 | | |

The two verifiers use python-flint. The first enumerates equal residues by combinations with replacement and tests $T$ itself. The second forms the initial sextuples from ordered halves of three, lifts each coordinate separately, and deduplicates the resulting multisets. It tests $6(T+1)$ directly and checks all prime claims with recursive Lucas primality proofs. The counts agree at every stage.

```sh
python3 check_seven_perfect_powers.py
python3 check_seven_perfect_powers_cross.py
```

## Power-residue conditions for an unrestricted partner

Use the [auxiliary primes and interval bound](factor_seven_residue_obstructions.md#uniform-counts-in-intervals-arbitrarily-far-from-the-origin). Given a large interval length $H$, take $h$ as the product of the first $J$ odd primes, with $J$ maximal subject to $7^h\le H$. Write $L$ for its largest prime. Then $L\sim\log\log H$. At every divisor $d>1$ of $h$ we have one distinct prime $r_d$ of exact order $d$ for 7.

**Outside at most $O_c(T_c(H))$ primes $p$ in any interval of $H$ integers, the conjecture holds at $1+7p^a q^b$ simultaneously for all odd primes $q$ and all $a,b\ge1$ for which some prime $\ell\le L$ divides $b$ but does not divide $a$.** Here $0<c<\log2$ and
$$T_c(H)=H\exp\!\left(-\exp\!\left(c\,\frac{\log\log H}{\log\log\log H}\right)\right).$$

For an odd $\ell\mid h$ select the $s=2^{J-1}$ divisors $d$ that are multiples of $\ell$. Then $\ell\mid r_d-1$. For $\ell=2$ select all $s=2^J-1$ divisors $d>1$. At each prime $r=r_d$, use the map $x\mapsto x^{(r-1)/\ell}$ from the units onto the $\ell$th roots of unity in the residue field. Every value has exactly $(r-1)/\ell$ preimages.

The folded exponents give at most $K=\binom{h+5}{6}$ possible residue vectors for $t$. Exclude the at most $2^J-1$ choices $p=r_d$ across all tests. If $\ell\mid b$ and $\ell\nmid a$, the power-residue value of $p$ at any nonzero coordinate is determined by the value for $t$ and the common choice of $a\bmod\ell$, because raising to $a$ permutes the $\ell$th roots of unity. At most one coordinate can be zero, necessarily $r_d=q$, giving at most $\ell$ free choices there. Thus at most $\ell(\ell-1)K$ vectors for $p$ can occur. Their product modulus divides $7^h-1\le H$. The Chinese remainder theorem bounds the number of integers in our interval per vector by $2H/\ell^s$.

Summing over $\ell$ bounds the exceptional count by
$$2^J+2H(J+1)L^2K\,2^{-2^{J-1}}.$$
Since $J\sim\log\log H/\log\log\log H$, this is $O_c(T_c(H))$. All estimates are uniform in the location of the interval and in the partner $q$.

If $b$ has a prime divisor $\ell\le47$, either $\ell\mid a$ and the perfect-power theorem applies, or the preceding interval estimate applies for all sufficiently large $H$. This proves the density-one statement in the introduction. The density in $b$ is $1-\prod_{\ell\le47}(1-1/\ell)$, with the product over primes.

**A growing range.** Fix $0<\epsilon<1$. For a set of primes $p$ of relative density one, every odd prime $q$, every positive $a$ all of whose prime divisors are at most 47 (including $a=1$), and every $b$ having a prime divisor at most $(1-\epsilon)\log\log p$ satisfy the conjecture at $1+7p^a q^b$. In particular this holds for all $2\le b\le(1-\epsilon)\log\log p$. To see this, apply the interval estimate on each dyadic interval $[H,2H)$, where $L\sim\log\log H$. A divisor $\ell\le47$ is already covered, and a larger one cannot divide the stated $a$. The number of exceptional $p$ in each interval is $o(H/\log H)$, so the prime number theorem gives the density claim.
