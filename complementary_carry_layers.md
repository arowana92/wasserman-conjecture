# Exact carry counts at two separated prime-power layers

Let $p$ be prime, $s,j\ge1$ and $d>1$ a divisor of $p^s+1$. For nonnegative integers $X_1,\ldots,X_m$ with total $X$, put
$$
 C_j=\lfloor X/p^j\rfloor-\sum_i\lfloor X_i/p^j\rfloor .
$$
This is the contribution of layer $p^j$ to the $p$-adic valuation of the multinomial coefficient. Choose signed integer representatives $u_i\equiv dX_i\pmod{p^{j+s}}$, and put $u=\sum_i u_i$. Suppose that every $|u_i|$ and $|u|$ is strictly less than $p^j$. Define
$$
 w_d(v)=\begin{cases}
 1& d\nmid v,\\
 0& d\mid v,\ v\ge0,\\
 2& d\mid v,\ v<0 .
 \end{cases}
$$
Then
$$
                     C_j+C_{j+s}=\sum_i w_d(u_i)-w_d(u). \tag{1}
$$
In particular, if the two parts of a binomial and their total have representatives nonzero modulo $d$, **exactly one of the two layers contributes a carry**. For three parts with all four representatives nonzero modulo $d$, the layer contributions sum to two.

## Proof

Set $M=p^j$, $B=p^s$, $h=(B+1)/d$. For any one part or total $Y$ with its representative $v$, write $dY=MBz+v$. We claim
$$
                \lfloor Y/M\rfloor+\lfloor Y/(MB)\rfloor=hz-w_d(v). \tag{2}
$$
Indeed the two quotients before taking floors are $Bz/d+v/(dM)$ and $z/d+v/(dMB)$. The perturbations have absolute value less than $1/d$. If $d\nmid v$, then $d\nmid z$ and $d\nmid Bz$. The two unperturbed fractional parts are nonzero multiples of $1/d$ and sum to one, since $B\equiv-1\bmod d$. Neither perturbed quotient crosses an integer. If $d\mid v$, coprimality of $d,MB$ gives $d\mid z$, so both unperturbed quotients are integers. Each floor decreases by one exactly when $v<0$. This proves (2). The parameters $z$ add from the parts to the total, proving (1).

The complementary fractional parts are the elementary identity behind Midy's theorem. Here the perturbed floors give a local valuation formula with signed errors. Lean checks the signed floor identity (2), including integral centers, and the resulting binomial identity when all three representatives are nonzero modulo the denominator. [Lewittes on Midy's theorem](https://arxiv.org/abs/math/0605182) · [Lean proof](ComplementaryCarryLayers.lean)

If the original congruences hold modulo $p^L$ with $L\ge j+s$, the same pairing holds at every pair of layer indices $i,i+s$ with $j\le i\le L-s$. In the nonintegral binomial case, every complete block of $2s$ consecutive indices in $[j,L]$ contributes exactly $s$ to the valuation. This statement counts the contribution of the specified layers and allows further carries elsewhere.

## Consequence for all three rows at a weight

Let $N=1+D p^k Q$ with $Q$ a power of a prime different from $p$, $k\ge1$, and $\gcd(D,pQ)=1$. Suppose $D$ is a prime power, twice an odd prime power, 12, or 15. The full conjecture holds if there exist $d>1$, $j,s\ge1$ and a signed integer $v$ such that
$$
 d\mid p^s+1,\quad dQ\equiv v\pmod{p^{j+s}},\quad D|v|<p^j,
 \quad d/\gcd(d,v)>D . \tag{3}
$$
The [two-power assignment theorem](two_prime_powers.md) forces a cut $aQ,(D-a)Q$ to add without carries in base $p$ in any failure, for some $0<a<D$. The representatives $av,(D-a)v,Dv$ satisfy (1), and none is divisible by $d$. Exactly one of the two indicated layers must carry. The same argument applies at any further cofactor where the two variable prime bases are forced onto the same row. [Derivation of the cut](paired_power_orders.md)

For any fixed odd prime $\ell>5$ congruent to 3 modulo 8, let $s=\operatorname{ord}_\ell(2)/2$. **For every $k\ge1$ the full conjecture holds at $1+5\cdot2^k\ell^f$ whenever $2^{s+1}\mid f+1$, $f\ge1$.** Since 2 is a quadratic nonresidue at $\ell$, its order is even and $2^s\equiv-1\bmod\ell$. For each even $n$ one has $v_2(\ell^n-1)=2+v_2(n)$. The odd part of the exponent leaves the valuation unchanged, $v_2(\ell^2-1)=3$, and each further squaring contributes one. Thus $\ell^{f+1}\equiv1\bmod2^{s+3}$. Apply (3) with $d=\ell$, $j=3$, $v=1$. All other odd prime bases in the factor-five family were already excluded by the [low-bit conditions](predecessor_five_seven.md#full-factor-five).

For example this covers $\ell=11$ at every $f\equiv63\bmod64$, and $\ell=43$ at every $f\equiv255\bmod256$, simultaneously for all binary exponents. In fact the forced binary cut carries specifically at layer $2^{s+3}$. At $2^3$ the residues of $2Q,3Q$ are $6,1$ since $Q\equiv3\bmod8$.

The same conclusion works for every fixed positive integer $b$ in place of 1 in the exponent shift. Take $d=\ell^b$, $s=\operatorname{ord}_{\ell^b}(2)/2$. The order is even and its half power is $-1$, the unique nonidentity square root of one modulo an odd prime power. The sufficient condition is now $2^{s+1}\mid f+b$. For odd $b$ this supplies a progression of odd exponents for each negative residue $-b$.

Run **python3 check_complementary_carry_layers.py** for signed representatives, vanishing residues, exact factorial comparisons and the sample prime-power families. Run **lean ComplementaryCarryLayers.lean** with Lean 4.36. [Arithmetic verifier](check_complementary_carry_layers.py)
