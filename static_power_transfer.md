# Four prime-power factors with two independently varying exponents

Let $a=p^i$ and $b=q^j$ be prime powers of positive exponent, $D=ab$, $n=ab+1$. Let $L_m$ denote the exponent of the unit group modulo $\operatorname{lcm}(1,\ldots,m)$. In every statement below $t$ is any positive integer with at most two distinct prime factors, and $e$ is **any positive multiple of $L_{n^3}$**. The primes in $t$ may also divide $abn$. Equivalently the last two prime powers may have independently chosen positive multiples of $L_{n^3}$ as exponents.

- **Every $a=2$ or $a=3$ works with every $b$.** Any three positive trinomials of weight $N=1+ab t^e$ share a prime.
- **Every $a$ works with every $q$ when $L_{pa}\mid j$.** The condition on $j$ depends only on the first prime power. The condition on the final two exponents depends only on the first two prime powers.
- At $a=4$, the assertion holds unless $b\equiv2\bmod3$ and $b+1,3b$ have disjoint binary ones. Among the 9,700 prime powers $b\le100{,}000$, these tests certify all but 11.
- For every fixed $a=p^i$, in any interval of $H$ integers at most $O_a(H^{\beta_p})$ prime powers $b$ can fail the assertion even once, where $\beta_p=\log_p(1/(2\sin(\pi/(4p+2))))<1$. The estimate also holds for $b=q^j$ with a fixed $q\ne p$ and $j$ in any interval of $H$ consecutive exponents, with a constant depending on $a,q$. At $p=2$, $\beta_2=0.69424\ldots$.

The proof reduces a putative failure to one of the following **binomial cuts at base $p$**, at the single reference weight $n$. At least one must add without a carry:
$$
 \begin{aligned}
 &cb,\ (a-c)b+1 &&(2\le c<a),\\
 &f n/w,\ (w-f)n/w &&(2\le w<a,\ w\mid n,\ 1\le f<w).
 \end{aligned} \tag{1}
$$
Thus the test requires $O(a^2)$ cuts for any $b$, independent of the size of the two remaining powers. [General reference-weight construction](two_power_seed.md)

## Both varying bases in one row

As in the reference construction, write $W=t^e$. The two companions to the endpoint have coprime scales $U,V$ multiplying to $N$. If both varying bases go to $A$ of scale $U$, choose a part $xW$ at the intersection of their full-power zero positions. Then $U\mid x$, $3\le U\le x<D$, and $U$ is a unitary divisor of $n$. This includes the cases $\omega(t)\le1$. Put $v=n/U\ge2$. At every prime of $D$, at least one of the small cut $(x,n-x)$ and the small companion $C_0=v(f_1,f_2,f_3)$ avoids that prime, with $\sum f_i=U$. In particular if the prime is a varying base, the cut avoids it. [Transfer at shared bases](two_power_seed.md#both-bases-in-the-same-row)

First suppose $p\ne q$. If $A$ takes $q$, then $x=mb+\delta$ with $1\le m<a$ and $\delta\in\{0,1\}$. Indeed if $q\mid W$, then $b\mid x$ by the full-power intersection. Otherwise $W\equiv1\bmod b$ and $xW$ has residue zero or one modulo $b$. The bound $3\le x<D$ gives the displayed range. The positive integer $U$ divides
$$
 a x-m n=a\delta-m,
$$
a nonzero integer of absolute value less than $a$. Consequently $U<a$. The cut $(x,n-x)$ has an unshifted part $cb$ with $1\le c<a$. Both cut entries are divisible by $U$, which is coprime to $b$, so $U\mid c$ and $c\ge3$. If the small companion takes $p$, use any of its cuts for the second line of (1) with $w=U$. Otherwise the first line works.

If $C$ takes $q$, then $q\nmid W$. Two of its normalized parts are multiples of $b$, giving $U\ge2b+1$. Hence $2\le v=n/U<a$. Since $V\equiv v\bmod b$, the same two entries of $C_0$ are positive multiples of $b$. Cut at one of these entries, $cb$. Here $c<a$ since two such entries and a positive third share the total $ab+1$, and $v\mid c$ since $\gcd(v,b)=1$. If $C_0$ avoids $p$, this gives the first line of (1). Otherwise $(x,n-x)$ avoids $p$, and $U\mid x$ gives the second line with $w=v$.

## Two different varying bases in different rows

Write $W=QR$ and the selected parts as $UQx,UQy$ and $VRz,VRh$. Positivity of the remaining parts gives $(x+y)(z+h)<D$. If a varying base divides $D$, absorb the rest of its power into $Q$ or $R$. The remaining static cofactor is a prime power or 1, which excludes a split by the four-factor test. Hence in a split with $p\ne q$ we can take $a,b,Q,R$ pairwise coprime. [Four-factor test](two_prime_powers.md)

Each of $a,b$ divides a displayed multiplier in its assigned row. They must go to the same row, as placement in different rows would give a product of sums at least $(a+1)(b+1)>D$. They divide different displayed multipliers there since a single multiplier is less than $D$. Relabel so $a\mid x,\ b\mid y$. At $p$ the $y$ part must have residue one modulo $a$ because it cannot be a multiple of $a$ as well as $b$. Thus $p\nmid y$. Likewise $q\nmid x$.

Put $g=\gcd(x,y)$, $X=x/g$, $Y=y/g$. We have $a\mid X,\ b\mid Y$. At each prime $p,q$ the ratio projection in the reference construction gives two positive entries $\ell X,\ell Y$ adding with a nonnegative third entry, carry-free at total $n$. Write their possibly different multipliers as $\ell_p,\ell_q$. The static congruences give $\ell_p Y\equiv1\bmod a$ and $\ell_q X\equiv1\bmod b$. Therefore
$$
 M=\ell_q X+\ell_p Y\equiv1\pmod{ab}.
$$
Both terms are positive, so $M\ge ab+1=n$. On the other hand $\ell_p,\ell_q\le n/(X+Y)$, so $M\le n$. Equality forces
$$
 \ell_p=\ell_q=k=\frac{n}{X+Y}.
$$
This step makes a single reference cut work at both static primes. It is verified in Lean along with the bounded-scale identity above. [Lean proof](StaticPowerTransfer.lean)

Since $g(X+Y)(z+h)<ab=n-1$, we have $k>g(z+h)\ge2$. Also $X+Y\ge a+b$, so $k\le(ab+1)/(a+b)<\min(a,b)$. The two reference parts are $kX,kY$ and sum to $n$. The second is $cb$ with $c\ge k\ge3$ and $c<a$, since $kX\ge3a>1$. They add without a carry at $p$, proving (1).

If $p=q$, the static cofactor $D$ is itself a prime power. After absorption of shared factors the four-factor test excludes a split. In a grouped assignment $U<D$ forces $A$ to take the base of $D$, since the companion would require two positive normalized parts divisible by $D$. If this base divides $W$, then $D\mid x$. Otherwise $W\equiv1\bmod D$ and $x\equiv0,1\bmod D$. Both contradict $3\le x<D$. This also proves the result for any single static prime power.

## Small first powers and a second specialization

For $a\le3$, a split is impossible because $3\le k<a$. In a grouped assignment $U<a$ is impossible, and so is $U\ge2b+1$ because $U$ is a proper divisor of $ab+1$ and $(ab+1)/2<2b+1$. This gives both unrestricted small cases.

For $a=4$, the same proper-divisor bound excludes $U\ge2b+1$. Thus the grouped case must have $U=3$, implying $b\equiv2\bmod3$. The congruence $3\mid4\delta-m$ forces $x=3b$ or $b+1$. The companion has three equal parts, so it carries at 2. The cut therefore specializes to the stated binary test. In a split $k=3$, and $b\mid Y$ with $Y<(4b+1)/3<2b$ gives $Y=b$. The very same test follows.

Now take $b=q^j$, $q\ne p$, $L_{pa}\mid j$. Then $b\equiv1\bmod pa$, and $pa=p^{i+1}>a+1$. The first line of (1) reduces in those low digits to $c,a-c+1$. Both are at least two, whereas the only two base-$p$ units in $a+1$ are $1,a$. For the second line, $w<a$ is coprime to $b$, so $b\equiv1\bmod w$ and $w\mid a+1$. Also $w$ is a unit modulo $pa$. The two entries now specialize modulo $pa$ to positive multiples of $(a+1)/w\ge2$ whose sum is $a+1$. Again the two units cannot be allocated. This proves the exponent theorem, including unrestricted choices of the two final prime bases.

## Counting the possible second powers

We give the digit estimate for the finite list (1). For fixed positive slopes and nonnegative integral offsets, the number of extendible residues modulo $p^m$ at which two affine functions of one nonnegative integer add carry-free is $O(\Lambda_p^m)$, with a constant depending on the slopes and $p$, where $\Lambda_p=1/(2\sin(\pi/(4p+2)))$. The bound holds with arbitrary input higher digits. Here are details for the offsets and the possibility that the sum slope has the largest valuation.

Use the digits of the two functions and the digitwise complement of their sum. Their three digits sum to $p-1$ in every column of an accepted addition. Remove the common $p$-power from the three slopes by shifting the words, rounding offsets down. The low columns are fixed. If they fail there is no accepted input. The minimum valuation now occurs at least twice. If all three slopes are units, varying one input digit after an accepted prefix permutes each of the three output digits. Their sum is $p-1$ or $2p-1$, since the lower columns have no incoming addition carry. Averaging gives $(p+1)/2$ accepted choices. This case requires odd $p$.

Otherwise one of the words has a delay $d>0$. After removing the delay, its digits give bijective input coordinates, including when this word is complemented. At position $j\ge d$, one of the undelayed output digits is at most $p-1$ minus the true coordinate at $j-d$. It is a permutation of the $p$ possible current input digits given the past. Encode the current coordinate by its zero-based rank among the allowed coordinates in increasing order. The rank is at most the true digit. Hence encoded digits a distance $d$ apart sum to at most $p-1$, and the encoding is injective. Each resulting chain of length $v$ has at most $p\Lambda_p^{v-1}$ strings by the operator norm of the matrix $(\mathbf1_{u+w\le p-1})_{0\le u,w<p}$. Its all-ones Rayleigh quotient is $(p+1)/2$. The first $d$ free positions cost only a fixed factor. [Matrix norm and digit encoding](inverse_digit_sieve.md#encoding-the-digits-of-the-cofactor) · [Binary affine bound](binary_affine_counts.md)

Apply this estimate to the first line of (1) directly. In the second line restrict $b$ to its required class $b_0\bmod w$, taking $0\le b_0<w$. In the parameter $(b-b_0)/w$ the two entries have positive integral slopes and nonnegative offsets. Taking the input modulus to be the first $p$-power at least $H$ bounds each count by $O_a(H^{\beta_p})$. There are only finitely many cuts for fixed $a$. If instead $b=q^j$ with $q$ fixed and $q\ne p$, its order modulo $p^m$ is $\gg_{p,q}p^m$. Thus a residue occurs only $O_{p,q}(1)$ times in an interval of $H$ exponents at that modulus. Division by $w$, a unit at $p$, preserves this bound on each allowed class. [Prime-power order bound](prime_power_ray_density.md#counting-exponents)

## Reproduction

Run `lean StaticPowerTransfer.lean` with Lean 4.36 and `python3 check_static_power_transfer.py`. Through $ab\le500$ the checker tests 6,190 primitive pairs, finds 80 exact cross-prime equalities, and cross-checks 219 exclusions with the general finite reference test. The digit estimate is checked against 600 exact prefix counts, including 184 with the sum slope delayed. It also checks the 9,700 prime powers $b\le100{,}000$. At $a=5,7,8,9$ it certifies respectively 9,602, 9,375, 9,627 and 9,565 of them for all admissible $t,e$. [Verifier](check_static_power_transfer.py)
