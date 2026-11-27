# Exponent classes for the 31-bit Mersenne prime

Set $R=2^{31}-1$. At $N=1+dR2^e$, for $1\le d\le512$, $e\ge1$, the following statements hold for **every** selection of three positive trinomial rows.

* All exponents work for 508 of the 512 multipliers. The four left open are $413,481,485,493$.
* Every multiplier through 512 works for all $e\le10^{30}$, and for every square or cube exponent without a size bound.
* At any one exponent, at most one multiplier in the entire range can fail. In any two consecutive exponents, at least one works for all 512 multipliers.

The [17 remaining arithmetic progressions](mersenne31_exponent_classes.txt) are necessary conditions for a failure. Each line gives $d,U,i,j,a,b,E,H$, with $e\equiv E\bmod H$. The search has five possible row shapes.

| $d$ | $U$ | $(a,b)$ | Exponent period $H$ | Classes |
|---:|---:|---:|---:|---:|
| 413 | 292096724993 | (35,59) | 393892433651712 | 4 |
| 481 | 344222663173 | (65,74) | 32012707674996 | 3 |
| 485 | 343614161929 | (50,97) | 13455207182160 | 4 |
| 493 | 284039356417 | (17,58) | 30818270171136 | 5 |
| 493 | 352724230657 | (58,68) | 37580345198760 | 1 |

## Enumerating the scales without bounding the exponent

Take $d$ odd first. Use the [two-row description for Mersenne primes](mersenne_multipliers.md). Write $W=2^e$, $UV=N$. The rows are
$$
 (UiW,UjW,(dR-U(i+j))W+1),\qquad V(Ra,Rb,U-R(a+b)).
$$
The scales $U,V$ consist of full prime powers of $N$ assigned to their respective avoiding rows. In particular $\gcd(U,V)=1$. The positive integers satisfy $i<j$, $a\le b$, $(i+j)(a+b)<d$, and $U$ is odd and coprime to $dR$. At each prime power $w\parallel d$ one row must have residues $0,0,1$ modulo $w$. For a trial $U\bmod d$, enumerate $W\in\langle2\rangle\bmod d$ and keep the residue if for every $w$ either $UiW,UjW$ are both zero or one is zero and the other one, or the two second-row residues $Ra/U,Rb/U$ satisfy the same test.

For each surviving shape, enumerate $R(a+b)<U<dR/(i+j)$ by branching on the bits of $U$ from low to high. The bits of $iU,jU$ computed so far must be disjoint submasks of $dR$. Reject an interval with no lift of the current prefix. Once at most one lift lies in the interval, check all its bits. In base $R$ the two parts $VRa,VRb$ must fit under the second digit of $N$, so
$$
              (aU^{-1}\bmod R)+(bU^{-1}\bmod R)\ \le\ d2^e\bmod R,
$$
using least residues. The right side depends only on $e\bmod31$. Exactly 42,252 choices of scale and shape reach a nonempty set of exponent residues. An independent enumeration through $d=127$ uses the unrestricted-height theorem and generates the binary scales with no adjacent ones, reproducing all 1,536 choices there.

At $p^f\parallel U$ strip $p^f$ from the first two parts when testing their sum. The lowest $f$ digits of the results are known since $W\equiv-(dR)^{-1}\bmod U$. Reject a carry. The condition $2^e\equiv-(dR)^{-1}\bmod U$ gives one class modulo $\operatorname{ord}_U(2)$ if soluble. Factor the order into prime powers $w$ and project to each order-$w$ subgroup by raising both sides to $\operatorname{ord}_U(2)/w$. The baby-step giant-step search covers every exponent from 0 through $w-1$. Recombine the logarithms and verify the resulting power modulo $U$. Compare this class by CRT with $e\bmod31$ and with every exponent modulo $\operatorname{ord}_d(2)$ that passes the residue and carry tests at the full powers of $d$. There are 73,940 joint classes. A separate logarithm calculation at the prime powers of $U$ gave the same subsequent list of 153 classes as the search using the order.

## Eliminating whole progressions by factorial valuations

In a class $e\equiv E\bmod H$, at a layer $w=p^j$ the possible $W$ form the orbit of $2^E$ under multiplication by $2^H$ modulo $w$. Its length is $\operatorname{ord}_w(2)/\gcd(H,\operatorname{ord}_w(2))$ for odd $p$. At $p\mid U$ the first row must avoid $p$. At $p\nmid U$ the second row can be reconstructed since $V\equiv(1+dRW)/U\bmod w$. If $p\mid N$, the second row is required to avoid $p$ by the assignment of full weight powers. If $p\mid N-1$, at least one row must avoid it. For any row with parts $z_i$, absence of a carry up to $w$ requires $\sum_i(z_i\bmod p^l)<p^l$ at every $1\le l\le j$.

The verifier first uses bounded orbits at the primes of $d$, small primes of $U$ and other primes through 97. If the orbit exceeds 128 elements it omits the test at that layer. This leaves 153 classes. For each required small prime put $H'=\operatorname{lcm}(H,\operatorname{ord}_p(2))$ and $w=p^{v_p(2^{H'}-1)}$. The possible initial residues $W\bmod w$ are supplied by $H'/H$ exponent representatives. Keep the ones passing the tests, then examine all $p$ lifts of each surviving residue at the next power. At $p^f\parallel U$ also require $p^{f+1}\nmid N$ once that layer is reached. Emptying a level excludes the class. Stopping after 35 levels or when the frontier exceeds 150 keeps the class as a possibility. The result is the table above, with $i=1,j=2$ in every case.

All four odd multipliers exceed 256, so no even multiplier through 512 can reduce to one of them. Each pair of distinct remaining classes is incompatible under CRT. The compatibility test also fails when one exponent is increased by one, including pairs from the same class. The classes omit squares and cubes, checked by exhibiting a divisor $v\le4096$ of the period on which the class residue is absent from the list of squares or cubes modulo $v$.

## The finite exponent check

In every remaining shape, all odd primes of $d$ are forced onto the second row because neither 1 nor 2 is divisible by them. Write $e=E+Ht$. The odd integer $d$ has exactly two prime factors $p,q$. For either prime $p$ we have $2^H\equiv1\bmod p$. The value $V=(1+dR\,2^{E+Ht})/U$ is fixed modulo
$$
       w_p p^j,\qquad w_p=p^{v_p(d)+v_p(2^H-1)}
$$
by $t\bmod p^j$. This follows from the binomial expansion of $(2^H)^{p^j}$, or induction on $j$. The algorithm lifts $t$ one digit at a time. It keeps only residues whose second row has no carry at all layers tested. It can stop at any completed level. An unfinished level is discarded in full, leaving the previous list as a valid necessary condition.

Let the resulting lists at the two primes be $S\bmod A$ and $T\bmod B$. To enumerate exactly the $t\le L=\lfloor(10^{30}-E)/H\rfloor$ that pass both tests, write $t=x+Ak$ for $x\in S$. Sort the residues $yA^{-1}\bmod B$ for $y\in T$. Then $k+xA^{-1}$ modulo $B$ belongs to this sorted list, with $0\le k\le(L-x)/A$. Binary search in the appropriate intervals enumerates precisely these $k$, including wraparound. Each candidate is checked by modular powering and the factorial-layer test on the second row at both primes, up to 80 base-$p$ digits. The respective counts excluded at this last stage for $d=413,481,485,493$ are **72,189, 103,756, 34, and 1,371,529**. There are no survivors. Direct modular checks also compare both individual prefix sieves for every $0\le t\le1024$ in every progression.

Use Python with python-flint and a C++17 compiler.

~~~sh
c++ -std=c++17 -O3 check_mersenne31_scales.cpp -o /tmp/wasserman-m31-scales
python3 check_mersenne31_exponents.py --scale-checker /tmp/wasserman-m31-scales --compare mersenne31_exponent_classes.txt
python3 check_mersenne31_range.py
~~~
