# Six digit units in the units place and two higher positions

Let $p\mid N-1$ be prime, with $s_p(N)=6$. Suppose the five positive-exponent digit units lie at just two positions $k,k+d$ with $k,d\ge1$. **The full conjecture holds if $p^d\le5{,}000{,}000$, for every $k$.** The result includes repeated units at either position. With just one positive position the [equal-gap theorem](equal_predecessor_gaps.md) already suffices.

Set $X=p^d$. Two disjoint pairs at a positive equal distance are necessary, so each position needs at least two units. The two pairs both have gap $d$. The remaining unit is at $k+\delta d$, where $\delta=0$ or 1. Its exponent must be an odd multiple of $d$. Thus only $k=nd$ with $n+\delta$ odd requires testing. Both positions carry at least two units and one carries three, so $p\ge5$. Write
$$
F=X+1,\quad m=2F+X^\delta,\quad N=1+mX^n,\qquad
A=(1+X^{n+\delta},X^nF,X^nF).
$$
A counterexample needs this row and a companion $C$ avoiding every prime of $m$ and also 2. Indeed $A$ carries at 2 because of its equal parts. Its assigned unitary weight divisor $U$ divides $F$, which in turn divides $N$. The normalized companion parts at scale $N/U$ sum to $U$.

## Testing all outer exponents together

For each odd full power $w=\ell^e\parallel F$, the two repeated entries divided by $w$ have the same residue
$$z=(-1)^n(F/w)\pmod w.$$
If adding this residue to itself produces a carry in the first $e$ base-$\ell$ digits, the row cannot avoid $\ell$. Omit this power and the entire power of 2 from $F$, calling the remaining product $H$. This gives $U\mid H$. Write $C=(N/H)(r_1,r_2,r_3)$, with positive $r_j$ of sum $H$. The scale is an integer, since $H\mid F\mid N$.

For each full power $w=\ell^e\parallel m$, the residues of the $r_j$ are $H,0,0$ in some order, because $N\equiv1\bmod w$. Thus $27m^2\le H^3$ by the product bound and $2w+1\le H$ by positivity. Since $H<m$, the Chinese remainder theorem determines each part uniquely after choosing its exceptional powers. Two further tests eliminate almost all of these rows.

- For any pair $a,b$ of normalized entries, $(a+b)/\gcd(a,b)$ divides the binomial coefficient of the corresponding scaled pair, since multiplying that coefficient by $a/(a+b)$ gives an integer. A common factor of this divisor and $2m$ excludes the row.
- At $w=\ell^e\parallel m$ the two entries $r_j$ divisible by $w$ have next blocks in the actual row given by $(r_j/w)H^{-1}\bmod w$. A carry between them at any of these $e$ digits excludes the entire family.

For the survivors, evaluate $X^n$ modulo 64 for every exponent of its required parity. At each residue the bottom six binary digits of the actual three entries are determined since $H$ is odd. At least one of those six columns carries in every case.

| Stage | Count |
|---|---:|
| Prime powers $X\le5{,}000{,}000$ of prime base at least five, counting the two choices of $\delta$ separately | 697,808 |
| Passing both size tests at $H$ | 103,856 |
| Unordered normalized companion rows | 616 |
| Excluded by a pair divisor | 432 |
| Further rows excluded by a predecessor digit block | 169 |
| Further rows excluded by the binary orbit | 15 |

The first calculation constructs single coordinates from subsets and matches their sums. The independent audit uses a sieve to factor integers, normalizes at the full odd part of $F$, partitions the predecessor powers among three slots, and checks factorial floors on the resulting actual entries. It recovers the same 616 rows after scaling and the same three exclusion counts. With python-flint installed, run both programs.

~~~sh
python3 check_two_level_units.py
python3 check_two_level_units_cross.py
~~~

## Large prime powers on the side with the lowest unpaired unit

The following restriction applies at arbitrary gaps and without the two-position hypothesis. In the [general equal-gap row](equal_predecessor_gaps.md), suppose $u\le i,j$, put $B=p^{i-u}+p^{j-u}$ and $F=p^d+1$. Then $m=(N-1)/p^u=BF+1$. **Every full prime power $w$ of either $F$ or $m$ must satisfy $w^2<BF$ in a counterexample.** In particular this holds at the two-position profile $\delta=0$ with $B=2$.

For $w=\ell^e\parallel m$, normalize the companion at scale $Z=N/F$. Then $Z\equiv F^{-1}\equiv-B\bmod w$. The two positive normalized entries divisible by $w$ are $wx,wy$, with $w(x+y)\le F-1$. Their next blocks in base $\ell$ are at least $w-Bx$ and $w-By$. Carry freedom implies their sum is at most $w-1$, so $w+1\le B(x+y)$. This gives $w^2<BF$.

For $w=\ell^e\parallel F$ assigned to $A$, divide its two unshifted entries by $w$. Since $p^u\equiv-1\bmod w$, their residues modulo $w$ are the negatives of $p^{i-u}F/w$ and $p^{j-u}F/w$. If $w^2\ge BF$ these two positive numbers sum to at most $w$. Their negative residues sum to at least $w$, forcing a carry. Finally if the base of $w$ is assigned to the companion, then $U\mid F/w$. The bound $27m^2\le U^3$ gives $w\le F/(3m^{2/3})<F^{1/3}$ and again $w^2<BF$.
