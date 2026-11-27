# Uniform exceptional counts on arithmetic progressions

Write $h(x)=-x\log_2 x-(1-x)\log_2(1-x)$, $a=h(1/6)$ and $b=h(1/3)$. Let $H$ be the number of terms of an increasing arithmetic progression in $[X,2X]$. In each of the following cases at most $O(HX^{-\delta})$ of its weights can fail Wasserman's full three-row conjecture, with some $\delta>0$ independent of the first term and the common difference.

1. **Every progression** with $H\ge X^\theta$ for fixed $\theta>b/(1+b-a)=0.72405194\ldots$.
2. **Every odd common difference $d$** with $dH\ge cX$ and $H\ge X^\theta$, for fixed $c>0$ and $\theta>a=0.65002242\ldots$.

There is no coprimality condition on the residue of the progression. For example, at each odd $d\le X^{0.3499}/2$ the conclusion holds uniformly in **every** class modulo $d$ throughout $[X,2X]$. The [previous interval estimate](two_predecessor_congruences.md#the-unrestricted-three-row-count) includes consecutive weights and certain prescribed binary digits.

## The arithmetic information in one part

In a counterexample choose a binary carry-free nonendpoint row $A$, assign 2 to it, and write $P=(U,S)$ for the full powers assigned to it from $N,N-1$; their complementary products are $V,T$. Put $K=US$. Every part $z$ of $A$ is at least three and satisfies $U\mid z$, $S\mid z(z-1)$. Thus $z$ allows at most $X^{o(1)}$ possible $P$, uniformly. Also $N$ is in one residue class modulo $K$. The complementary normalized row has sum $U$ and part product divisible by $T^2$. Hence
$$ K\ge3(N-1)^{2/3}S^{1/3}. $$
The number of weights in the block admitting any such assignment with $K\le Y$ is $O(YX^{o(1)})$ uniformly for $1\le Y\le X$. These facts follow from the [endpoint theorem](all_weights_endpoint_reduction.md) and the [two-congruence reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product).

Let $F=\gcd(K,d)$. In $H$ progression terms a specified $P$ permits at most $1+HF/K$ weights. On the other hand each $z$ is zero modulo $\gcd(U,d)$ and zero or one independently at each prime power exactly dividing $\gcd(S,d)$. We may stratify by both gcds and all three residues modulo $F$ at total cost $X^{o(1)}$, by the divisor bound applied to $d$. For fixed residues, deleting the bottom $f=\lfloor\log_2 F\rfloor$ bits of any part loses no information: the remaining prefix places it in $2^f\le F$ consecutive integers.

## The two portions of the progression that vary

Set $t=v_2(d)$, $u=\lceil\log_2 d\rceil$, $l=\lfloor\log_2 H\rfloor$, and pad the binary weights to length $L=\max(\lceil\log_2(2X+1)\rceil,u+l)$, so $L=\log_2 X+O(1)$. The bottom $t$ bits are fixed, and the bits above position $u+l-1$ have $O(1)$ possible joint values. Fix those values. In the intervening positions use these two words:

- The $u-t$ bits starting at position $t$ run in an odd-step progression modulo $2^{u-t}$. Each word occurs $O(H/2^{u-t})$ times since $2^{u-t}\ll X/H\ll H$ in both ranges.
- The $l$ bits starting at position $u$ also have bounded multiplicity. The integers $\lfloor N/2^u\rfloor$ progress by zero or one, each at most twice, and their range has length at most $H$. Reduction modulo $2^l$ costs only a constant factor.

Fix $\epsilon>0$. Except at $O(HX^{-c_\epsilon})$ terms, each of these words and each of its suffixes obtained by deleting low bits has at most half its length plus $\epsilon L$ ones. Under a uniform word the probability of this upper tail on $m>0$ bits is at most $\exp(-2\epsilon^2 L^2/m)$ by the binomial tail bound. The two multiplicity bounds transfer it to the progression, and a union over $O(L)$ suffixes preserves a power saving.

## Counting the remaining exceptions

Choose canonical witnesses and put the uniform distribution on the remaining exceptional weights in one nonempty stratum. Use base-two Shannon entropy $\mathcal H$. Since $N$ determines $P,A$, $A$ determines $N$, and $\mathcal H(P\mid A_i)=o(L)$, conditional subadditivity gives
$$
 \mathcal H(N)\le\tfrac13\sum_i\mathcal H(A_i)+\tfrac23\mathcal H(N\mid P)+o(L).
$$
Indeed the sum of the three conditional part entropies given $P$ is at least $\mathcal H(N\mid P)$. This is the [three-part entropy inequality](common_factor_count.md#information-that-every-individual-part-supplies).

On the remaining fixed positions after deleting $f$ bits, a known one can belong to only one of the three parts, costing at most $b$ bits on average by entropy concavity. A known zero costs zero. On the remaining variable positions the total expected number of ones across all three parts is at most half the number of positions plus $2\epsilon L$, since the parts partition the ones of $N$. Concavity and the tangent to $h$ at $1/6$ bound the mean cost here by $a$ times the number of positions plus $O(\epsilon L)$. Thus if $K\ge K_0$ throughout the stratum,
$$
 \mathcal H(N)\le
 b\big(L-u-l+(t-f)_+\big)
 +a\big(l+u-\max(t,f)\big)
 +\tfrac23(l+f-\log_2 K_0)_+ +O(\epsilon L)+o(L).
 \tag{1}
$$

For the first assertion discard $K<H X^{-\eta}$ using the small-$K$ count, at cost $O(H X^{-\eta+o(1)})$. The last term involving $K_0$ in (1) is then at most $2(f+\eta L)/3+O(1)$. The resulting expression is convex in $0\le f\le u$ because $b>a$. At $f=0$ its main part is at most $b(L-l)+al$. At $f=u$ it is $b(L-l-u)+al+2u/3$, which has the same upper bound since $b>2/3$. For the stated $\theta$, this upper bound is below $l$ by a positive multiple of $L$. Take $\eta,\epsilon$ small enough and sum over the strata.

In the second case $t=0$ and $u+l=L+O_c(1)$. Use the unconditional $K_0\gg X^{2/3}$ directly. Apart from errors the upper bound is now $a(L-f)+\tfrac23(l+f-2L/3)_+$. By convexity it suffices to check $f=0,u$. At zero it is $aL+\tfrac23(l-2L/3)_+$, strictly below $l$ with a linear margin when $l/L\ge\theta+o(1)>a$ (use $a<2/3$). At $u$ it is $al+2L/9+O_c(1)$, likewise below $l$ since $a(1-a)>2/9$. This completes both proofs.

Run `python3 check_progression_weights.py` for independent integer checks of the progression blocks, the reconstruction after deleting low bits, and the numerical entropy margins.
