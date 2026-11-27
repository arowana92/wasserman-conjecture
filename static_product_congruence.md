# A polynomial modulus for products of ordinary prime factors

Fix distinct prime powers $a=p^i$, $b=q^j$ with positive exponents. Use $D=ab$, $n=ab+1$, $c=b^{-1}\bmod a$, $d=a^{-1}\bmod b$, and $g=\gcd(c,d)$, taking least positive inverses. For $s=p,q$ let $B_s$ be the least power of $s$ strictly greater than $Dn$. Set
$$
 M_0=B_p B_q,\qquad M=\operatorname{lcm}(B_p,B_q,c,d).
$$

**Theorem.** Suppose at least one of $g\ge3$, $p\nmid\binom b d$, $q\nmid\binom a c$ fails. For every positive integer $W$ with at most two prime divisors and $W\equiv1\bmod M$, the three-row conjecture holds at $1+DW$. In particular $W$ may be the product of two ordinary primes whose individual residues are arbitrary inverses of each other modulo $M$. This gives an infinite progression of primes for the second factor for every fixed first prime coprime to $M$, by Dirichlet's theorem. The modulus satisfies $M<D^4(D+1)^2$ since $B_s\le s Dn$, $pq\le D$ and $cd<D$.

If either binomial test fails, the smaller modulus $M_0$ suffices. More precisely, under just $W\equiv1\bmod M_0$, a failure requires both binomial tests and at least one of
$$
 g\ge3,\qquad \gcd(1+DW,c)\ge3,\qquad \gcd(1+DW,d)\ge3 . \tag{1}
$$
Compare the [canonical cut for two independently chosen large exponents](canonical_cofactor_cut.md). That result also permits shared bases. The present modulus is a condition on the product and allows exponent one on both factors.

## Projection of a grouped assignment

Write $N=1+DW$. The [endpoint assignment](two_prime_powers.md) gives two nonendpoint rows of coprime scales $U,V\ge3$ with $UV=N$. Call the row taking the varying bases $A$ if those bases go to one row. This includes $W=1$ or a single prime power. Some part of $A$ is $xW$ with $3\le U\le x<D$ and $U\mid x$. The other row is $V(f_1,f_2,f_3)$ with $\sum f_i=U$.

If $A$ takes a static prime $s$, the cut $xW,N-xW$ is carry-free there. Reduction modulo $B_s$ gives the actual reference entries $x,n-x$, since $W\equiv1$ and their sum is $n<B_s$.

If the other row takes $s$, let $r_i$ be the residues of $V f_i$ modulo $B_s$. They sum to $n$ without carries, and
$$
 U r_i\equiv n f_i\bmod B_s .
$$
Both products are at most $Dn<B_s$, hence equal. Write $h=\gcd(f_1,f_2,f_3)$, $u=U/h$, $f'_i=f_i/h$. Then $u r_i=n f'_i$ for all $i$, so $u\mid n$ by primitivity of the $f'_i$. The positive reference row $C_0=(n/u)(f'_1,f'_2,f'_3)$ is the same at every static prime taken by this row. Its normalized sum $u$ divides $x$ and is at least three. Lean checks the lift from the congruence to equality and the primitive-denominator step. [Lean proof](CanonicalCofactorCut.lean)

If this row took both static primes, some entry of $C_0$ would be divisible by $ab$ in a positive triple summing to $ab+1$, an impossibility. If it took just $p$, then $A$ would take $q$. At $b$ we have $x=mb+\delta$, $1\le m<a$, $\delta=0,1$. Hence $u\mid ax-mn=a\delta-m$ and $u<a$. But $C_0$ avoiding $p$ requires two positive $f'_i$ divisible by $a$ since $n/u$ is coprime to $a$. This is impossible. Exchanging $p,q$ handles the other mixed choice.

Therefore $A$ takes both static primes. Its cut $x,n-x$ equals $bc,ad$ up to order by the congruences at $a,b$. The cut avoids both primes, and $U\mid x$ with $\gcd(U,ab)=1$ forces $U\mid c$ or $U\mid d$. This proves the grouped alternatives in (1), without projecting $U$ when the companion takes no static prime.

## Projection of a split assignment

Write $W=QR$ with the two full prime powers assigned to different rows. Both are coprime to $D$. The [four-factor argument](two_prime_powers.md) places both static powers on one row, in distinct multipliers of $UQx,UQy$, with $a\mid x$, $b\mid y$. Indeed placement on different rows would give a product of multiplier sums at least $(a+1)(b+1)>D$, and placement on the same multiplier would put it at least at $D$. Write $s_0=\gcd(x,y)$ and $X=x/s_0,Y=y/s_0$. Since $UQy\equiv1\bmod a$ and $UQx\equiv1\bmod b$, we still have $a\mid X$, $b\mid Y$.

At each static layer the residues of $UQx,UQy$ are positive since $UQ$ is a unit there and $x,y<D$. They fit into $n$ because $N\equiv n\bmod B_p,B_q$. Their cross products by $x,y$ are at most $Dn$ and congruent, hence equal. At the two primes this yields reference pairs $l_p(X,Y), l_q(X,Y)$ with sums at most $n$, and $l_pY\equiv1\bmod a$, $l_q X\equiv1\bmod b$. Thus $l_qX+l_pY$ is congruent to one modulo $ab$, lies above one and is at most $n$: it equals $n$. The weighted average of $l_p,l_q$ then equals their common upper bound $n/(X+Y)$, so both equal $k=n/(X+Y)$. The reference cut is $kX=ad$, $kY=bc$. It avoids both static primes. Finally $(x+y)(z+t)<D$ for the two positive multipliers $z,t$ on the other row, so $g=k>s_0(z+t)\ge2$. This proves (1).

When $W\equiv1\bmod M$, we have $N\equiv n$ modulo both $c,d$. As $n=bc+ad$, $\gcd(c,n)=\gcd(d,n)=g$: use $\gcd(c,a)=\gcd(d,b)=1$. All three alternatives in (1) now reduce to $g\ge3$.

## A polynomial period, including shared prime bases

The same canonical test certifies every $N=1+Dt^e$, $\omega(t)\le2$, for **every positive multiple $e$ of $\lambda(M)$**, where $\lambda(M)$ is the exponent of the units modulo $M$. In particular the sufficient period is less than $D^4(D+1)^2$. Each prime in $t$ can independently have any positive multiple of this period as exponent.

To see this, call the full varying power of each base in $t^e$ its varying factor. At any divisor $v$ of $M$, each varying factor coprime to $v$ is one modulo $v$. For a split with neither base shared with $D$, this gives all the static congruences in the preceding proof. A split with a shared base is impossible: after absorbing the shared static power into its varying factor, the four-factor criterion is applied with at most one remaining static prime power.

In a grouped assignment, a shared base $s$ is taken by $A$. Its full power in $D$ divides $x$ because the part $xW$ was chosen at the intersection of the zero positions at the *full* predecessor powers of the varying bases. Write $W=s^v Z$ with $s\nmid Z$. Above the shift $s^v$, the cut at $xW$ gives the carry-free sum $xZ+(D-x)Z$ in base $s$. Since $Z\equiv1\bmod B_s$, this projects to $x,D-x$. Both units digits are zero, so adding one to the second entry gives the required reference cut at $n$. All other static primes use the preceding projection unchanged. Thus both static primes go to $A$ and $U$ divides $c$ or $d$. Since $U\mid M$ and $U\mid N$, it is coprime to $W$. Every varying factor is one modulo $U$, giving $U\mid n$ and then $U\mid g$. This proves the period claim also when the varying bases divide $cd$.

## Reproducible example

Take $a=5,b=7$. Then $c=d=3$ and the base-5 test fails since $\binom73=35$. It suffices that $W\equiv1\bmod 5^5 7^4=7{,}503{,}125$. For example $103\cdot582767=8\cdot7{,}503{,}125+1$ gives $N=2{,}100{,}875{,}036$. A separate certificate for this weight uses the four primes $7,563,71761,2100875033$. No row avoids two of them, so one divides all three coefficients. The verifier enumerates the 9,072 sorted rows avoiding 7 and the two avoiding $N-3$ and checks the other valuations by factorial floors. A row avoiding both 563 and 71761 must have all parts divisible by their product, leaving only the positive triples of total 52 to inspect. Run **python3 check_canonical_cofactor_cut.py** for this certificate and the exact low-layer specialization of the companion on 14,202 finite instances, including 4,904 where the original scale does not divide the seed. Run **lean CanonicalCofactorCut.lean** with Lean 4.36 for the two projection lemmas. [Verifier](check_canonical_cofactor_cut.py)
