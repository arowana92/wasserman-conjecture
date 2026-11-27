# Almost every prime works up to its tenth root in both parameters

**For a set of primes $q$ of relative natural density one, Wasserman's full three-row conjecture holds at every**
$$
 1+dq^e,\qquad 1\le d,e\le q^{1/10}.
$$
The number of exceptional primes in $[X,2X]$ is $O(X^{0.997})$. The parameters are simultaneous and their bounds grow with $q$.

There are also estimates uniform in an arbitrarily large fixed exponent or base. Put $\alpha=\log_2\lambda=0.551463\ldots$, where $\lambda^3=\lambda^2+1$, $\lambda>1$. Let $1\le D\le H$.

- In any increasing arithmetic progression of $H$ positive integers, at most $O(D^{4-\alpha}H^\alpha)$ prime terms $q$ can fail at some $1+dq^e$, $d\le D$, for a given $e\ge1$. The constant is absolute, independent of $e$, the first term, and the step.
- For every fixed odd prime $q$, in any increasing arithmetic progression of $H$ positive exponents, at most $O(D^{4-\alpha}H^\alpha)$ terms can fail for some $d\le D$. The same uniformity includes the prime $q$.
- For the 1,000 cofactors $1\le d\le1000$ the corresponding bounds in both directions are $O(H^{0.336})$. At $d=986$ alone both improve to $O(H^{0.284})$. In particular, for almost every prime $q$, all 1,000 cofactors work for all $e\le q^{0.66}$. At 986 every $e\le q^{0.71}$ works for almost every prime. These statements also apply to every cubic exponent in the latter range for all $d\le1000$, by the [cube reduction](denominator_residue_hulls.md).

## Reading the digits that actually vary

Let a carry graph in base $p$ permit at most $C p^{s k}$ possible next $k$ digits from any state reached by a common low prefix, when the word must extend to a finite accepting word. Consider a sequence of integer inputs on indices satisfying
$$
 v_p(Y_i-Y_j)=h+v_p(i-j)\quad(i\ne j),\qquad h\ge0.
$$
They have a common prefix modulo $p^h$. Among $M$ consecutive indices the next $k=\lceil\log_p M\rceil$ digits distinguish every index. Thus at most $C p^s M^s$ indices can pass. The bound does not depend on the possibly very large $h$.

Two applications give the needed uniformity. For $Y_i=(A x_i^e+B)/U$ with $p\nmid U$, $x_i$ on an arithmetic progression of positive units, split indices modulo $p$ (modulo 4 at 2). More precisely retain the unit classes after splitting the original progression this way. On each class the input step is $S$ divisible by $p$ (by 4 at 2). For any distinct inputs in the class,
$$
 v_p(x_i^e-x_j^e)=v_p(x_i-x_j)+v_p(e).
$$
For $Y_i=(Aq^{a+Si}+B)/U$, with $q>1$ a unit, instead split indices modulo $p-1$ (modulo 2 at 2). Then the new step $S'$ makes $q^{S'}$ congruent to 1 modulo $p$ (modulo 4 at 2) and
$$
 v_p(q^{a'+S'i}-q^{a'+S'j})=v_p(q^{S'}-1)+v_p(i-j).
$$
In both formulas multiplication by $A$ adds $v_p(A)$ and the unit denominator preserves the valuation wherever the quotient is integral. The identities follow by factoring a difference of powers for exponents prime to $p$ and applying the binomial expansion for raising a principal unit to the $p$th power. At 2 the congruence modulo 4 ensures the exact one-step lift. These arguments still distinguish actual integral candidates if other indices give only $p$-adic integers. Each splitting uses at most a fixed number of classes depending only on $p$.

In the binary three-summand case with positive slopes $s_1,s_2,s_3$, include $s_0=s_1+s_2+s_3$ and put $J=\max_i v_2(s_i)-\min_i v_2(s_i)$. From any compatible low prefix, the remaining multiplications have the same slopes and new nonnegative offsets. The [binary affine digit theorem](binary_affine_counts.md) gives the preceding bound with $s=\alpha$ and $C=\lambda^J$, independent of the offsets.

## Summing over growing cofactors

For $Q=q^e$ a failing weight has the [endpoint row](all_weights_endpoint_reduction.md). A row avoiding $q$ has the form $A=(aQ,bQ,cQ+1)$, $a+b+c=d$. Assign full prime powers of $N$ to avoiding rows, with product $U$ assigned to $A$ and product $V$ to the other row $C$. Then $C=V(f,g,h)$ where $U=f+g+h$, $UV=1+dQ$. Each of $A$'s parts is divisible by $U\ge3$. Hence $a,b,c\ge1$, $U\mid a,b$ and $U\le d$. At least one of $A,C$ must avoid 2. When $A$ carries at 2, $U$ is odd: if $N$ is even its full 2-power must be assigned to $C$. [Assignment details](binary_affine_counts.md#consequences-for-wassermans-conjecture)

For the weights $w(f,g,h)=\lambda^J$ in the binary bound,
$$
 \sum_{f+g+h\le T} w(f,g,h)\ll T^3
$$
over positive triples. Indeed $w\le\sum_{t\in(f,g,h,f+g+h)}2^{\alpha v_2(t)}$, counting the four entries with multiplicity. Each entry takes any given value $O(T^2)$ times and $\sum_{t\le T}2^{\alpha v_2(t)}\ll T$, by a geometric series since $\alpha<1$. Thus for binary-free $A$ the total cost over $d\le D$ is $O(D^3 H^\alpha)$ by the preceding index bound.

For binary-free $C$, fix odd $U\le D$. In the progression of candidate bases split the indices by the residue of $q$ modulo $U$. There are $T\le U$ cyclic classes, where $T$ is $U$ divided by the gcd of $U$ with the progression step. Discard nonunits. In a progression of exponents at fixed $q$, instead use the period $T\le U$ of $q^e$ modulo $U$; a divisor $U$ of $N$ is coprime to $q$. In either case the fixed residue of $Q$ forces one residue $d\equiv-Q^{-1}\bmod U$. This leaves at most $1+D/U\le2D/U$ choices of $d$ per class. There are at most $2H/T$ indices per class. At these indices apply the binary prefix argument to $V=(dQ+1)/U$ and the slopes $(f,g,h)$. For each triple the cost is
$$
 \ll T(D/U)(H/T)^\alpha w(f,g,h)\le D U^{-\alpha} H^\alpha w(f,g,h).
$$
Group the triples by dyadic ranges of their sum $U$. The weighted triple bound above gives $\sum_{U\le D}U^{-\alpha}\sum_{f+g+h=U}w(f,g,h)\ll D^{3-\alpha}$. Adding the two rows' bounds gives the claimed $O(D^{4-\alpha}H^\alpha)$. A possible prime term $q=2$ in the base progression adds at most one.

For primes in $[X,2X]$ take $H\asymp X$, $D=\lfloor(2X)^{1/10}\rfloor$ and sum the estimate over all $e\le(2X)^{1/10}$. The resulting exponent is $1/2+9\alpha/10<0.997$. Indeed $\lambda<733/500$ and $(733/500)^{125}<2^{69}$. Division by the prime number theorem estimate for the number of primes gives relative density one. More generally, for $d\le q^u$ and $e\le q^v$ with $0\le u\le1$, the same conclusion holds whenever $\alpha+(4-\alpha)u+v<1$.

## The 1,000 bounded cofactors and ratios of exceptional exponents

In the [bounded enumeration](prime_base_density.md#counting-exponents-for-any-fixed-prime) one of the forced carry graphs admits a positive vector $w$ with $Bw\le p^{0.336}w$, where $B$ is its adjacency matrix. This bounds paths from every state by $O(p^{0.336k})$. There are finitely many graphs. Their tested primes are distinct from the prime base of a failing weight: the small prime bases are already settled. Applying the digit-shift argument in both directions gives the two uniform bounds asserted above. Summing over exponents gives the density claims since $0.336+0.66<1$. At the two rows of denominator 986 the productive base-seven graphs admit exponents $s=0.284,0.248$ respectively. The verifier constructs positive integer vectors and checks the inequalities raised to the power 1000, giving the asserted $0.284$ bound. [The two possible rows](cofactor_986_powers.md)

These last exponents also bound rational ratios *simultaneously across different prime bases*. In one of the two rows, every exceptional $W$ satisfies $v_7(W-1)=1$. Modulo $7^{k+1}$ write $W=8^z$, where $z\bmod7^k$ is a unit and 8 has order $7^k$. There are $O(7^{sk})$ eligible $z$ by the path bound. If $W_1,W_2$ are powers of the same prime giving failures, their reduced exponent ratio $a/b$ must satisfy $a/b\equiv z_1/z_2\bmod7^k$ by raising to the powers $b,a$. The exponents are units at 7. There are at most $O(7^{2sk})$ ratio residues. In particular the number of positive integer multipliers $r$ in any interval of $H$ integers for which **any** pair $1+986q^e,1+986q^{re}$ could both fail is $O(H^{0.568})$, by taking $7^k\asymp H$ and summing over rows. The [single-row theorem](cofactor_986_powers.md) applies to every pair.

For rational ratios of height $\max(a,b)\le T$, choose $T^2<7^k\le7T^2$. Distinct reduced positive ratios of this height have different residues, by the cross determinant. The counts of ratios allowed by this necessary test are therefore $O(T^{1.136})$ in the first row and $O(T^{0.992})$ in the second. Combined with the exclusion of the first row at 29th powers this gives an infinite consequence:

**There is an effectively specified set $\mathcal P$ of relative density one among the primes such that, for every prime $q$ and positive $E$, at most one prime $j\in\mathcal P$ can give a failure among all $1+dq^{87Ej}$, $d\le1000$.** To construct it, for each pair of distinct prime multipliers $a>b$ discard $a$ if the second-row ratio test passes. Use $T$ the smallest power of two at least $a$ to select the precision. Discard 7 as well, and test only units at 7. The number discarded through $X$ is $O(X^{0.992})$ on summing the preceding height bound over dyadic ranges. Distinct prime multipliers already give reduced ratios, so no failing pair can survive this construction. At these exponents only 986 remains by the cube reduction, and [the 29th-power obstruction](common_power_exponents.md) excludes its first row. The set $\mathcal P$ is independent of $q,E$.

Run the exact digit and arithmetic controls with Python and python-flint.

~~~~sh
python3 check_binary_affine_counts.py
python3 check_moving_power_bases.py
~~~~
