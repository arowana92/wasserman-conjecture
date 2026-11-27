# All exponents at bounded odd binary cofactors

Write $N=1+d2^a$, with $d$ odd and $a\ge1$. For $1\le d\le65535$ the full three-row conjecture holds at **every exponent** for 32,758 of the 32,768 cofactors. It holds for **all** the cofactors whenever $a\le10^{24}$. In particular every odd $d<12031$ works without an exponent bound.

The following are the only progressions left by the all-exponent proof. In every interval of $K\ge2$ consecutive exponents, at most $O(K^{0.36})$ exponents can give an exception at **any** of the 32,768 cofactors. The last column gives the exponent in the corresponding individual upper bound. In the table $q=2^a$, the row avoiding 2 is $(xq,yq,zq+1)$, and $U,T$ are defined below.

| $d$ | $(x,y,z)$ | $U$ | $T$ | $a$ modulo $H$ | $H$ | carry base | count exponent |
|---:|---|---:|---:|---:|---:|---:|---:|
| 12031 | 681, 9222, 2128 | 3 | 1 | 907 | 5876 | 3 | .14 |
| 15871 | 2065, 13450, 356 | 5 | 1 | 6002 | 7772 | 5 | .23 |
| 33661 | 820, 32840, 1 | 5 | 1 | 410 | 820 | 5 | .26 |
| 56279 | 4676, 16513, 35090 | 7 | 1 | 801 | 1743 | 7 | .29 |
| 56279 | 4676, 33026, 18577 | 7 | 1 | 717 | 1743 | 7 | .29 |
| 60829 | 2065, 57736, 1028 | 7 | 1 | 35049 | 89610 | 7 | .28 |
| 61423 | 16962, 34177, 10284 | 11 | 1 | 7030 | 9520 | 11 | .36 |
| 62463 | 8712, 53603, 148 | 121 | 47 | 5401 | 24310 | 3 | .12 |
| 62809 | 20544, 42264, 1 | 3 | 1 | 15529 | 31058 | 3 | .11 |
| 65095 | 31301, 32770, 1024 | 113 | 47 | 312 | 644 | 5 | .22 |
| 65243 | 16536, 36930, 11777 | 3 | 1 | 14788 | 31980 | 3 | .13 |

## From binary digits to two prime-power assignments

At a counterexample the endpoint row $E=(1,1,N-2)$ occurs, with coefficient $N(N-1)$. Exactly one of the two other rows need be chosen to avoid 2. Call it $A$. Its binary digits are disjoint, so its parts have the stated form with disjoint $x,y,z$ summing to $d$. The two unshifted slopes can be ordered $0<x<y$. Assign each full prime power of $N$ to a nonendpoint row avoiding its prime, giving $UV=N$ with $U$ assigned to $A$. Then $U\mid x,y$ since $q$ is a power of two, and the companion is $C=V(f,g,h)$ with $f+g+h=U\ge3$. The shifted part is divisible by $U$, so $z>0$. We have $\gcd(U,d)=1$. If $d$ has at most one distinct prime factor, the existing two-prime predecessor theorem applies. [Endpoint reduction and the two-prime theorem](all_weights_endpoint_reduction.md)

Write $d=ST$ as a product of full prime powers, assigned to $A$ and $C$ respectively. If $w=p^j\parallel S$, the actual parts of $A$ are $0,0,1$ in some order modulo $w$. At least one of $x,y$ is zero modulo $w$. If both are zero, there is no extra restriction. Otherwise the other slope $k$ is a unit and $kq\equiv1\pmod w$. At $w\parallel T$ a full prime power, the residues of $(f,g,h)$ are $0,0,U$ in some order, since $V$ is a unit and $UV\equiv1\pmod w$. Therefore
$$2w+1\le U,\qquad T^2\mid fgh,\qquad 27T^2\le U^3.$$
Also $S\mid xy$, hence $d/\gcd(d,xy)\mid T$. These inequalities are only used to discard impossible assignments.

Two independent programs obtain exactly the same **658** tuples $(d,x,y,z,U,T,a_0,H)$. The last two entries describe the solution $a\equiv a_0\pmod H$ of all required power congruences, including $q\equiv-d^{-1}\pmod U$. The C++ program partitions each $d$ by submasks. It tests companion existence by choosing an exceptional coordinate at every prime of $T$ and combining the residues by CRT. The least positive representatives fit precisely when their sum is at most $U$, since the remaining room is a multiple of $T$. The Python program starts with $x+y$, splits its bits, and fills $z$ from the unused positions. It constructs individual companion coordinates and checks that their exceptional-prime sets are disjoint and that the three coordinates sum to $U$. It also computes the power congruences at the individual prime powers of $U$.

## Higher carries and interval counts

For each prime dividing $US$, the row $A$ must avoid it. At $p\mid U$ require $p^{v_p(U)+1}\nmid N$ as well. For each prime dividing $T$, the same companion must work, and $C$ must avoid that prime. Test every $q$ in the orbit of $2^{a_0}$ under multiplication by $2^H$ modulo the largest power of $p$ through 65536 (at least $p$ and at least $p^{v_p(U)+1}$ when $p\mid U$). At each lower power $w$ require the factorial floor of the total divided by $w$ to equal the sum of the three corresponding floors. There are 17 surviving tuples.

For each take the least prime $p\mid US$. Read the digits of $q$ from low to high. From the multiplication carries $c=(0,0,1)$ initially, an input digit $b$ forms the three numbers $c_i+b s_i$, where $s=(x,y,z)$. Their residues must sum to less than $p$, and their quotients by $p$ form the next state. There are finitely many states because each coordinate stays between 0 and $s_i$. Mark all states whose carries already add without overflow, then mark all their predecessors recursively. These are exactly the states admitting a finite completion. Six tuples have no permitted first digit leading to a marked state, leaving the eleven table entries.

For a surviving tuple let $P$ be the transition matrix on marked states reachable after an admissible first input digit, counting each distinct digit. Put $v=P^{32}\mathbf1$. Every component is positive since a successful completion can be padded with zeros. For each table entry the verifier proves componentwise, using integer powers, that
$$P^{32}v\le p^{32\delta}v,$$
where $\delta$ is the displayed count exponent. Iteration bounds the number of possible length-$j$ low prefixes by $O(p^{j\delta})$. Put $b=\operatorname{ord}_p(2)$ and $t=v_p(2^b-1)$. At $j\ge t$, the order modulo $p^j$ is $b p^{j-t}$, by the binomial expansion or lifting the exponent. Choosing $j=\log_p K+O_p(1)$ so the order is at least $K$ makes distinct exponents in an interval of length $K$ give distinct prefixes. Summing the eleven bounds proves the claim uniformly in the interval's location.

For the finite exponent bound, refine each progression by the multiplicative order of its step $2^H$ at each next $p$-power, following only marked carry edges. Each child fixes one additional digit of $q$. Stop when the period exceeds $10^{24}$, at which point each class contains at most one positive exponent in range. This produces 6,244,206 terminal exponents counted across the tuples. For every one, some forced prime has a carry in its first 80 digits, checked with modular powers and factorial floors. All computations use integers.

## The first cofactor where one odd row covers every lowest predecessor layer

Among the 658 tuples, the least $d$ with $T=1$ is 1903. This minimum concerns assignments with $U\ge3$ and a positive companion. For $d$ with at most one odd prime factor, a single row cannot cover the whole predecessor either. With at most two full powers in $N-1$, some part would miss the exceptional residue at every power and would be divisible by $N-1$, exceeding its maximum $N-2$.

The minimum is attained. For $a=119+2580j$, $j\ge0$, put $q=2^a$, $N=1903q+1$, $A=(264q,519q,1120q+1)$ and $C=(N/3,N/3,N/3)$. Here $3\parallel N$. The row $A$ has odd coefficient and its residues at each of $q,11,173$ are $0,0,1$ in some order. The full powers of $N$ split between $3$ on $A$ and $N/3$ on $C$. Thus all assigned lowest layers work. At the very next ternary layer the residues of $A$ modulo 9 are $6,3,3$. The prime 3 divides all three coefficients. At $a=119$ its exact valuation on $A$ is 81.

## Reproduction

From the repository root, with Python 3.9+ and a C++17 compiler:

```sh
mkdir -p .work
c++ -std=c++17 -O3 check_binary_cofactor_templates.cpp -o .work/binary_templates
.work/binary_templates > .work/binary_templates.txt
python3 check_binary_cofactor_rays.py --compare .work/binary_templates.txt
```

The Python verifier also runs independently with no options. Its optional `--from-templates` switch replays the carry and exponent calculations from a saved enumeration.
