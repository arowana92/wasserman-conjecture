# A factorial progression for arbitrary cofactors

Fix integers $B\ge2$, $J\ge0$ and $M\ge(J+3)B^2$. For **every prime** $p$, integer $u\ge1$ and cofactor $1\le d\le p^u\le B$, take any positive exponent
$$
                      e\equiv-u\pmod{M!}.
$$
Wasserman's full trinomial conjecture holds at $N=1+dp^e$. In fact, for every multinomial coefficient $K$ of weight $N$ with at least three positive parts,
$$
           p\nmid K\quad\Longrightarrow\quad
           N\left(\prod_{\substack{q\le B\text{ prime}\\ q\ne p}}q\right)^J\mid K. \tag{1}
$$
Thus a single factorial gives simultaneous progressions for every eligible base, a whole range of cofactors and any prescribed small-prime multiplicity. For example set $B=p^{2r}$ and $M=3B^2$. For every positive multiple $T$ of $M!$, **each of the $r+1$ consecutive exponents $T-2r,\ldots,T-r$ works for all $1\le d\le p^r$**, for every $r\ge1$.

### The valuation at a small prime

Put $b=p^u$ and $Q=p^e$. If $p\nmid K$, addition of its parts is carry-free in base $p$. The low $e$ digits of $N$ represent 1, so all but one part are multiples of $Q$. Take any two such parts $xQ,yQ$. They satisfy $x,y\ge1$ and $x+y\le d\le b$.

Fix $q\le B$ prime, $q\ne p$. Set $E=\lfloor\log_q(2B)\rfloor$, $l=\varphi(b)$ and $K_0=l(J+E+1)$. Then $\varphi(q^{K_0})\mid M!$. Indeed $qE\le q^E\le2B$, since $E\ge1$. Thus $qK_0\le B((J+1)B+2B)\le M$, so $M!$ has at least $K_0$ factors of $q$ counted with multiplicity. Also $q-1\mid M!$ and is coprime to $q$. We obtain $bQ\equiv1\pmod{q^{K_0}}$ by Euler's theorem. Consequently
$$
                    bN\equiv b+d\pmod{q^{K_0}}.
$$
Here $q^{K_0}>2B$ since $K_0\ge E+1$, so $v_q(N)=v_q(b+d)\le E$.

For each $j=1,\ldots,J+E+1$ use the layer $w=q^{jl}$. We have $b\mid w-1$ and $k=(w-1)/b\ge1$. Since $b(-k)\equiv1\pmod w$, we have $Q\equiv-k\pmod w$. The residues of the two parts are $w-kx,w-ky$, whose sum exceeds $w$ because $k(x+y)\le w-1$. Thus Legendre's floor formula gives at least one to $v_q(K)$ at **each distinct layer**. In particular $v_q(K)\ge J+v_q(N)$. The [Lean lemma](InverseScaleCarry.lean) checks the residue inequality over the integers.

### Large primes and the endpoint

If $q>B$ divides $N$, it divides neither $Q$ nor $x$, since $x\le d\le B$. The identity $(xQ)K=N K'$, where $K'$ is the multinomial after decreasing that part by one, gives $v_q(K)\ge v_q(N)$. Together with $p\nmid N$ this proves (1).

By the [endpoint reduction](all_weights_endpoint_reduction.md) it suffices to consider the endpoint with two other rows. If both coefficients are divisible by $p$, this prime divides all three. Otherwise (1) makes one of them divisible by $N$. For any remaining positive row with parts $c_i$, its coefficient is divisible by $N/\gcd(c_1,c_2,c_3)$ by the same decrement identity at every prime. This integer is at least three and divides the endpoint coefficient $N(N-1)$.

The construction also works inside any prescribed arithmetic progression of exponents at a fixed $p,d$. Choose $u$ large enough in the negative of the desired residue class, with $p^u\ge d$, take $B=p^u$, and increase $M$ to at least the requested modulus. Then every positive $e\equiv-u\pmod{M!}$ belongs to that progression. To place a finite block along a progression, choose its negative starting exponent far enough to leave $p^u\ge d$ at the last term.

Run `python3 check_inverse_scale.py` to verify the congruences and floor contributions, compare the multinomial assertion with exact coefficients at tractable exponents and check the [six small-offset progressions](close_power_gaps.md#a-short-progression-at-every-inner-exponent). The verifier also checks 15,036 layer inequalities and 18 certificates at inner exponents 257, 512 and 1001. With Lean 4.36.0 run `lean InverseScaleCarry.lean`.
