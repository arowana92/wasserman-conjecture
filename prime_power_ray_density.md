# Density one along each fixed odd-cofactor prime-power ray

**Theorem.** Fix a positive integer $d$ and a prime $q$ coprime to $d$. Among exponents $1\le k\le K$, at most $O_{d,q}(K^{1-\delta_d})$ can admit a counterexample containing $E=(1,1,N-2)$ at $N=1+dq^k$, for some $\delta_d>0$. If $d$ is odd, the bound holds for **all** counterexamples. In particular the full conjecture holds at a set of exponents of natural density one on every ray $1+d2^k$, including even $d$ after absorbing its power of two into the exponent.

The full assertion for odd $d$ follows from the [prime-power cofactor theorem](prime_power_cofactors.md): eventually any counterexample must contain $E$. An exceptional-set estimate among all integer weights alone gives no such estimate on these exponential sequences.

## A base-$p$ obstruction

Put $Q=q^k$ and suppose $E$ is in a counterexample. Its coefficient $N(N-1)$ is divisible by $q$. By Kummer a row avoiding $q$ has parts $(aQ,bQ,cQ+1)$ with $a+b+c=d$ and $a,b>0$. Assign all full powers of $N$ to avoiding rows. The product $U$ assigned to this row is at least 3, since the third row would otherwise have three positive normalized parts of sum at most 2. Thus some prime $p\mid U$ divides $a,b$, with $p\le d$ and $p\ne q$. The sum of $aQ,bQ,cQ+1$ must have no carry in this base.

Take $L=1+\lceil\log_p d\rceil$. If $Q$ has a block of $L$ consecutive digits equal to $p-1$ in base $p$, both $aQ$ and $bQ$ have digit $p-1$ at the last (most significant) position of that block. In fact the incoming multiplication carry for $aQ$ is an integer $j\in[0,a-1]$. After $L-1$ maximal digits the carry is
$$\left\lfloor\frac{a(p^{L-1}-1)+j}{p^{L-1}}\right\rfloor=a-1.$$
The following digit is $(a(p-1)+a-1)\bmod p=p-1$. The same holds for $b$. The two digits already force a carry. Hence at least one prime $p\le d$, $p\ne q$, must have its specified block missing from $q^k$.

## Counting exponents

For fixed $p\ne q$ the multiplicative order of $q$ modulo $p^l$ is $\gg_{p,q}p^l$. To see this, let $e$ be the order modulo $p$ for odd $p$, or modulo 4 for $p=2$, and put $v=v_p(q^e-1)$. For every positive $b$,
$$v_p(q^{eb}-1)=v+v_p(b).$$
The binomial expansion gives an increase of exactly one for each raising to the $p$-th power (using $v\ge2$ at $p=2$), and no increase for an exponent coprime to $p$. Thus the order modulo $p^l$ is $e p^{l-v}$ for $l\ge v$.

Choose $l$ as the largest multiple of $L$ with $p^l\le K$. Each residue of $q^k\bmod p^l$ occurs only $O_{d,p,q}(1)$ times for $k\le K$ by the order bound. Partition the $l$ low digits into blocks of length $L$. A string missing the forbidden block has at most $(p^L-1)^{l/L}$ choices, even if only these aligned blocks are tested. Therefore the number of exponents for this prime is
$$O_{d,p,q}(K^{\theta_p}),\qquad \theta_p=\frac{\log(p^L-1)}{L\log p}<1.$$
Sum over the finite set of possible $p$ and take $1-\delta_d$ to be the maximum of the displayed exponents over primes through $d$ (allowing $q$ in this maximum). If no such prime is available, the endpoint is already impossible. This proves the bound.
