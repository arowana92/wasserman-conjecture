# Two consecutive weights near a product of powers of one prime

**Theorem.** Fix any prime $p$ and integers $j\ge a\ge1$. For every $1\le d\le1000$, Wasserman's three-row conjecture holds at **both**
$$p^a(p^j+d)\quad\text{and}\quad p^a(p^j+d)+1.$$
The offset may be divisible by $p$. The powers and the two weights are unbounded.

Write $P=p^a$, $H=p^j$, $Q=H+d$, and for now let $N=PQ+1$. The following necessary condition for a counterexample holds even for $a>j$ and $d\ge H$:
$$27 Q^2<d^5. \tag{1}$$

## Proof of the bound

By the [all-weight endpoint reduction](all_weights_endpoint_reduction.md), a counterexample consists of the endpoint and two rows $A,C$. Assign each prime dividing $N(N-1)$ to a row among $A,C$ whose coefficient it does not divide. Call the row receiving $p$ the first row $A$. Write $UV=N$ and $PBT=N-1$ for the products of the full powers assigned to $A,C$, respectively $U,V$ and $PB,T$. Thus $Q=BT$, $\gcd(U,V)=\gcd(PB,T)=1$, and $U$ is coprime to $PQ$. When $p\mid Q$, its full power in $Q$ is included in $B$.

Kummer's criterion says that the row assigned a prime adds without carries in that base. Consequently $U$ divides all parts of $A$. The row $C$ divided by $V$ has three positive parts of sum $U$, so $U\ge3$. At every full power $w$ in $T$, two of these normalized parts are divisible by $w$. Thus each $2w+1\le U$, and the arithmetic–geometric mean inequality gives
$$27T^2\le U^3. \tag{2}$$

The parts of $A$ are $P\ell_i+\delta_i$, where exactly one $\delta_i$ is 1. Each $\ell_i$ is positive because each part is a positive multiple of $U$. The $\ell_i$ add without carries in base $p$ to $Q$. Since $Q\ge H$, some $\ell_i$ receives a digit unit worth at least $H$. Designate that part large. The remaining normalized parts $x,y$ satisfy $t=x+y\le d$. Choose $y$ unshifted, so $U\mid y$ and $U\le y<d$.

Factor $B=ef$ by putting in $e$ those powers whose residue 1 occurs at the actual part $Py$. The power of $p$ in $Q$, when present, goes into $f$ and divides $y$, since the entire $p$-power in $PQ$ divides $Py$. Every other power of $f$ also divides $y$, and $\gcd(f,U)=1$. Hence $f\mid y/U$ and $e\mid Py-1$. If the extra unit is on the large part then $e\mid x$. If the unit is on $x$, the large part is unshifted and $e\mid Q-t$. As $e\mid Q$, this gives $e\mid t$. In both cases $e\le d$, so
$$Q=efT\le \frac{dy}{U}\sqrt{U^3/27}<\frac{d^{5/2}}{\sqrt{27}}.$$
This proves (1), and $U\ge3$ gives $d\ge4$.

## Exhausting the offsets up to 1,000

For $j\ge a$ and $d\ge H$, both weights are at most $2d^2+1$, inside the [initial interval through one billion](full_range_1b.md). In the remaining cases we can restrict the search to $d<H$ and weights above that interval.

When $j=a$, $P=Q-d$, so $e\mid dy+1$. Enumerate $d$, $3\le y<d$, $3\le U\mid y$, $e\mid dy+1$ with $e\le d$, and $f\mid y/U$. For $B=ef$ the permissible $Q$ are units modulo $U$ satisfying $Q+Q^{-1}\equiv d\pmod U$, since $N\equiv0\pmod U$. Thus one can enumerate $T\equiv Q B^{-1}\pmod U$ subject to (2), then recover the actual integers $Q=BT$ and $P=Q-d$. Keep prime powers $P$ and check the unitary conditions. If the large part is shifted then $Ue\mid x$. Otherwise $t$ is a multiple of $e$ and $t\equiv Q\pmod U$. These progressions exhaust both orientations.

For $j>a$, (1) bounds $H$ directly and $P<H$. Factor $N$ and $Q$ for each weight. Enumerate all possible unitary $U<d$ and $T$ meeting (2) and the individual full-power bounds. It is necessary that $27Q^2\le Ud^4$ and $BU<d^2$, since $BU\le dy<d^2$. For each multiple $y$ of $U$ below $d$, recover $e=\gcd(B,Py-1)$ and use the same progressions for the other part. In both branches test every first row for carries at $p$ and every base dividing $UB$.

- When $j=a$, there are 119,644 assignments after the factor and size filters, producing 2,754 first rows. None passes all its assigned primes.
- When $j>a$, there are 109,392 weights above the initial interval before the assignment tests. There are 30 assignments and two first rows. Neither passes.

For the preceding weight $N_0=PQ$, the [prime-to-offset case](close_prime_powers.md) already covers $p\nmid d$. For $p\mid d<H$, put $q=p^{v_p(d)}$. In the normalized data $P'=Pq$, $H'=H/q$, $d'=d/q$, the offset is now coprime to $p$ and $d'<H'$. The size bound in that note, valid regardless of the order of the two exponents, yields
$$432(N_0-1)^2q\le P^3d^4.$$
Since $N_0-1\ge PH$ and $H\ge P$, this gives $P\le d^4/432$ and $H^2\le Pd^4/432$. Enumerate all powers satisfying these bounds with $p\mid d$, then apply the earlier per-weight factor-assignment and carry test to the normalized data. This leaves 28,383 weights, 6,534 assignments, and 146 first rows passing. None admits a complementary row even at the required prime-power residues. This proves the paired theorem.

## Reproduction and independent checks

With Python 3 and python-flint installed, run

```sh
python3 check_successor_powers.py
python3 check_successor_powers_cross.py
```

The second program reverses the square-case enumeration. For each prime power $P$ and each $d\le60$ it directly factors $Q$, extracts eligible full powers of $N$ by trial division, and tests all assignments above weight 10,000. Its set of 119 tuples $(P,d,y,U,B,e)$ equals the set obtained by solving the progressions. A separate check computes the trinomial coefficients directly with integer binomials for $P\le42$. All 65,256 first rows avoiding $p$ satisfy the derived divisibilities and appear in the predicted part progressions, including 15 with $U\ge3$ and 84 with $e>1$. The latter check takes the maximal assignments available to the first row and exercises the progressions even when the complementary size tests would reject them.
