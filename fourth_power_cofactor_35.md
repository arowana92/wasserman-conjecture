# Fourth powers at cofactor 35

Let $t$ be a positive integer with at most two distinct prime factors, and write $t=3^c s$, $3\nmid s$. The full three-row conjecture holds at $N=1+35t^4$ in each of the following ranges.

- $c\le5$, with no restriction on $s$ apart from the stated hypotheses.
- The sum of the ternary digits of $s$ is at most 16, with no size bound on $s$ or $c$. Through digit sum 17 the only remaining $s$ is the prime $21\,523\,361$, and it is excluded through $c=10^{12}$.
- $s<3^{32}=1\,853\,020\,188\,851\,841$ and $c\le10^{12}$. In this entire $s$ range, only the primes $21\,523\,361$ and $926\,510\,094\,425\,921$ still require an upper bound on $c$.

A second result concerns all offsets $|d|\le10000$. If $2s=3^n+d$, $n\ge1$, then the only $s$ that can occur at a failure have the form $(3^{2^j}+1)/2$, $j\ge4$. There is a further restriction on $c$ for these inputs below.

The [cofactor-35 reduction](cofactor_35.md) already handles even $t$, $c=0$ and prime-power $t$. In the remaining case $s$ is an odd prime power, $\gcd(s,105)=1$, and $t\equiv\pm3$ modulo both 7 and 8. Put $Q=3^{4c}$ and $R=s^4$. There is an odd-coefficient row $A$ of weight scale four that avoids all four prime bases of $35QR$. Its part $16QR$ is one modulo both 5 and 7, and is zero at $Q,R$. The other two parts are multiples of 35. Their exceptional residues at $Q$ and $R$ must be on different parts, since otherwise a positive part would be divisible by $35QR=N-1$. Thus
$$
 A=(16QR,\ 35zR,\ (19Q-35z)R+1),\quad
 0<35z<19Q,\quad 4\mid z,\quad 35zR\equiv1\bmod Q.        \tag{1}
$$
For the strict upper bound, the third part is a positive multiple of four. At 3, grouping the last two parts and dropping the first $4c$ digits gives the necessary carry-free addition
$$
                           16s^4+19s^4=35s^4.               \tag{2}
$$

## An unrestricted digit-sum theorem

Among **all** positive $s$ prime to 3 with ternary digit sum at most 13, exactly $3281$ and $265721$ pass (2). Both have two different prime divisors, since they equal $17\cdot193$ and $41\cdot6481$. Neither is eligible for the remaining case of $t$.

To see how to certify an unbounded input range, read $s$ one trit at a time, beginning with a unit. At each step discard a prefix $v$ if its digit sum exceeds 13 or its two residues $16v^4,19v^4$ at the new modulus add to at least that modulus. Once the budget is used, the prefix is the whole input and can be tested exactly. There is one prefix that also needs a proof for arbitrarily many intervening zero trits before the next nonzero digit.

That prefix is $b=3281$, of digit sum nine. Suppose $s=b+3^h u$, $h\ge40$, $u>0$ a unit of digit sum at most four. Since $35b^4<3^h$, the next $h$ columns in (2) require the addition of $16C u,19C u$ to have no carry through those columns, where $C=4b^3$. This follows from expanding the fourth power modulo $3^{2h}$. A separate four-unit prefix tree needs just the first six trits to force $u=92$. It checks a fully spent digit budget at its actual integer value and tracks the first failing column, which is also within that bound for every rejected terminal. For $u=92$, all five coefficients $35\binom4j b^{4-j}u^j$ are below $3^{40}$. Hence the five blocks of the expansion are separated and the top one requires (2) at $s=92$. It carries in the second ternary column. This rules out every such extension of $b$.

With this one stopping rule the whole prefix search visits 389,454 nodes and has no live prefixes at depth 70. Exact tests of exhausted budgets detect carries through position 80. A second search starts at each unit and jumps directly to the next occupied position, using the first failing modulus of the current integer as its bound. Both return exactly the two stated inputs.

To extend the theorem to 17 units, consider both ends of the ternary expansion. With a budget of eight units the low-prefix tree has no solution and forces a carry by the 39th place, including the first failures of all budget-exhausted leaves. Thus any success of length at least 39 spends at least nine units in its first 39 places. Indeed, if there were at most eight, truncation there would give a positive unit of that budget that passes all 39 places.

The leading 42 places likewise require nine units. Starting with each nonzero leading trit, use the real interval test for $[h,h+1)$ given below and keep every extension with total digit sum at most eight. The tree tests 14,385 prefixes and leaves no admissible 42-trit prefix. A length of at least 81 therefore requires 18 units because the two windows are disjoint. Through length 80 a direct ternary prefix search with budget 17 visits 21,703,037 nodes, stopping at exhausted budgets or length 80 and testing the remaining integer in full. It returns exactly $3281,265721,21523361$. The last one, $(3^{16}+1)/2$, has fifteen trits equal to one followed by a units trit equal to two, giving digit sum 17. This proves the stronger digit-sum conclusions.

## A finite range in the 3-free part of the input

For positive odd $s<3^{32}$ prime to 3, (2) holds at exactly seven integers:
$$
                       s=(3^n+1)/2,\qquad n=8,12,16,20,24,28,32 .
$$
There are two independent necessary prefix conditions. For the low trits of $s$, extend a residue modulo $3^k$ by one digit and test the two fourth-power residues at $3^{k+1}$. There are exactly two children at every clear unit prefix. In fact the new digit of $s^4$ ranges over all three choices by $3\nmid4s^3$. The output digits for $16s^4,19s^4$ and their sum reduced modulo three each permute all three digits. The differences of the digit totals over the three choices sum to three, so one choice carries.

For a high prefix $h$ with $j$ unused trits, $s/3^j\in[h,h+1)$. Choose a power of three $F\ge19((h+1)^4-h^4)$. If (2) holds then the quotients of $16s^4,19s^4$ by $F3^{4j}$ must themselves add freely. For multiplier $a=16,19$ its quotient lies between $\lfloor ah^4/F\rfloor$ and $\lfloor(a(h+1)^4-1)/F\rfloor$. Each interval has at most two members. Discard $h$ if no pair passes, and recursively extend the other prefixes. Taking potentially incompatible choices of the two quotients only enlarges this filter.

Using 19 leading and 23 trailing trits gives 1,839,223 high prefixes and 8,388,608 low residues. Intersect their intervals and residue classes and retain odd inputs. This yields 130,636,566 candidates, of which 672,164 pass the first 36 columns. Testing every column gives the seven values above. The program computes the exact fourth power in base $3^{24}$ to avoid large-integer overflow. Repeating with 18 leading and 24 trailing trits gives the same seven values. Direct arbitrary-precision enumeration agrees through $3^{13}$.

At $n=8,12,20,24,28$ the numbers have a proper divisor respectively $17,41,41,17,41$ coprime to its complementary factor. The two other numbers are prime by full-order certificates, using witnesses 6 and 7 respectively and the complete factorization of each input minus one.

For either prime $s$, the restrictions on $t\bmod 7,8$ prescribe $c=c_0+6i$ in one class with $1\le c_0\le6$. Write $W=81^c s^4$. The cut $16W,19W+1$ must add without a carry at 2, 5 and 7. A prefix of $i$ of length $k$ bits fixes $W\bmod2^{k+5}$ since $v_2(81^6-1)=5$. Start at $W=81^{c_0}s^4$ and branch, multiplying by $(81^6)^{2^k}$ for the occupied bit. Prune when $16W$ and $19W+1$ share any of the known low bits. Stop when the index is the only representative of its prefix with $c\le10^{12}$. At $p=5,7$ the scale modulo $p^7$ is periodic in $i$ with period dividing $p^6$, since $81^6=1\bmod p$. Test each of these seven layers at both primes and finally the low 110 bits of the cut.

| $s$ | Terminal indices | Passing the tests at 5 and 7 | Passing 110 binary places |
|---:|---:|---:|---:|
| 21,523,361 | 5,383,615 | 17,548 | 0 |
| 926,510,094,425,921 | 5,384,016 | 17,647 | 0 |

## Holding the 3-exponent fixed instead

For $1\le c\le5$, grow all residues $r=R\bmod Q$, starting at $1\bmod3$ and retaining every prefix that passes the cut $16r,19r$. Set $z=(35r)^{-1}\bmod Q$ in $[1,Q)$. Keep the size and divisibility conditions from (1). For each resulting shape, every possible $R$ is $1\bmod16$, an odd fourth power. Test the three affine parts through 12 binary places for all lifts of that prefix. For the shapes that survive, check every possible finite binary continuation by an exact graph.

Write $(X,Y,Z)=(16Q,35z,19Q-35z)$. From the multiplication carries $(u,v,w)$, an input bit $e$ produces $u+Xe,v+Ye,w+Ze$. The transition exists if at most one is odd, and the new carries are their integer halves. Seed the graph with the carries after the four low bits, including the $+1$ in the third part. A finite completion exists exactly when some reachable state has pairwise disjoint remaining carries, because supplying zero bits from that point reads those carries themselves as the remaining outputs. Enumerate distinct reachable states until the set closes. None accepts.

| $c$ | Candidate $z$ | Shapes needing the graph after 12 places | Largest graph |
|---:|---:|---:|---:|
| 1 | 2 | 0 | 0 |
| 2 | 26 | 2 | 19 |
| 3 | 301 | 23 | 3,131 |
| 4 | 4,501 | 413 | 7,955 |
| 5 | 71,002 | 6,748 | 10,190 |

## Why one pattern keeps appearing

For odd units $s=(3^n+d)/2>0$ with $|d|\le10000$, apply the [rational polynomial valuation formula](rational_polynomial_carries.md) to the polynomials $16(X+d)^4,19(X+d)^4$ and their sum, all divided by 16. The period is four. Choose $m$ with $3^m>1+\max_i35\binom4i|d|^{4-i}$. Test $1\le n<m$ directly and the first two eligible exponents on each residue class at $n\ge m$. Over all offsets the only zero valuation that persists is $d=1$, $n=0\bmod4$. There are no isolated zeros on the other stabilized classes. The small exponents yield seven additional representations of the same value $3281$. Thus as values, every survivor in these offset families equals $(3^{4k}+1)/2$, $k\ge2$.

There is also a direct check that this unbounded sequence really passes (2). For $X=3^{4k}$ with $k\ge2$, the base-$X$ blocks of $16s^4$ from most to least significant are $(1,4,6,4,1)$. Those of $19s^4$ are
$$
 1,\ (3X+61)/16,\ (15X+113)/16,\ (X+63)/16,\ (13X+19)/16.
$$
They are nonnegative and below $X$. Their residues modulo nine are $1,1,2,0,4$ respectively. The first list uses only the bottom two trits of each block, and every displayed pair adds freely there.

If $(3^n+1)/2$ is a prime power for $4\mid n$, then $n$ is a power of two. Otherwise put $n=2^u h$ with odd $h>1$ and $x=3^{2^u}$. An odd prime $\ell$ dividing $(x+1)/2$ would be the sole prime of $(x^h+1)/2$. For odd $h$, binomial expansion (or the elementary lifting-the-exponent lemma) gives $v_\ell((x^h+1)/(x+1))=v_\ell(h)$. But $(x^h+1)/(x+1)>h$ for $x\ge81$, $h\ge3$. The input at $n=8$ has two prime factors, leaving only $n=2^j$, $j\ge4$ in the offset theorem.

For these remaining prime powers put $4c=kn+b$, $0\le b<n$. The conjecture also holds whenever $k=0$, $k$ is odd, or $b=0$ and $k=2$. If $4c\le n$ then $16R=(3^n+1)^4\equiv1\bmod Q$. Consequently $35z=LQ+16$ in (1), with $0\le L<19$ (here $Q\ge81$). But $Q\bmod35\in\{1,11,16\}$, giving $L\bmod35\in\{19,24,34\}$, impossible. For odd $k$, $Q\equiv-3^b\bmod s$, since $3^n=2s-1$. As $4\mid b$ and $b\le n-4$, $35\cdot3^b<s$. The required cut $16Q,19Q$ at the base prime of $s$ already carries across the block of size $s$ because its residues sum to $2s-35\cdot3^b>s$. The cut is obtained from (1) by grouping and removing the shift $R$. Finally if $4c=2n$, $Q=(2s-1)^2$. The residues of $16Q,19Q$ at $s^2$ are $s^2-64s+16$ and $s^2-76s+19$, producing a carry since $s>140$.

## Reproduction

From the repository root run

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_fourth_power_cuts.cpp -o .work/check_fourth_power_cuts
.work/check_fourth_power_cuts 32 19 23
.work/check_fourth_power_cuts 32 18 24
python3 check_sparse_fourth_power_cuts.py
python3 check_fourth_power_scales.py
~~~

The sparse and offset checks use arbitrary-precision Python integers. The scale checker independently reverses the ternary inverse enumeration through $c=3$ and recomputes every terminal that passes the two static primes by direct modular powering and factorial quotients.
