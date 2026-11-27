# Prescribing the small-prime part of a common trinomial divisor

Let $m$ be any fixed positive integer and $B\ge3$ a cutoff including its prime factors. **At $\gg_{m,B}X^{0.55146}$ distinct weights through $X$, there are three positive trinomial rows for which the part of the coefficient gcd supported on primes through $B$ is exactly $m$.** More precisely, if $G$ is the gcd of their coefficients, then
$$
                 G>1,\qquad \prod_{p\le B}p^{v_p(G)}=m .
$$
One row is the endpoint $(1,1,N-2)$. The other two pass the lowest-layer avoidance tests at **every** prime dividing $N(N-1)$. For odd $m$, one of these two coefficients is odd. In particular any prescribed odd prime can be the sole common prime up to any prescribed cutoff, with exact valuation one. For $m=1$ all common primes exceed the cutoff.

The exponent in the count can be taken to be $\delta=\log_2\lambda$, where $\lambda^3=\lambda^2+1$ and $\lambda>1$. Even for a growing cutoff $B=\lfloor c_m\log\log X\rfloor$, the number is at least $X^\delta\exp(-O_m(\sqrt{\log X}))$ for some $c_m>0$. The construction specifies particular triples. The general three-row conjecture remains open.

## Keeping all small primes out of a primitive sum

First take $m$ odd. Let $a$ be the largest power of three at most $B$, and put $z=2^a$, $u=1+z+z^2$. No prime $p\le B$ divides $u$. The assertion is clear at 2, and at 3 we have $z\equiv-1$ because $a$ is odd. At any other odd $p\le B$, the order of $2$ modulo $p$ divides $p-1$ and all its factors of three are removed on raising to $a$. Thus $z$ cannot have order three. If $z=1\bmod p$, then $u=3\bmod p$, also nonzero.

Choose a positive odd $V$ whose binary ones have no two positions differing by $a$ or $2a$. Then
$$
              N=uV,\qquad A=(V,zV,z^2V)
$$
gives an odd trinomial coefficient. The three shifted bit sets are disjoint. We will impose $V\equiv1\pmod u$.

Let $W=2u\prod_{3\le p\le B}p$, with the product over primes. Choose relatively prime odd integers $d,e>u$ coprime to $W$, for example $d=W+1$ and $e=Wd+1$. Put $D=de$ and fix a unit $r_0$ modulo $D$ by
$$
                  ue r_0\equiv1\pmod d,\qquad ud r_0\equiv1\pmod e.
$$
Require $uV\equiv1+2D r_0\pmod{D^2}$. Then $r=(uV-1)/(2D)$ is integral and $r\equiv r_0\pmod D$, since $uV$ is odd. In particular $\gcd(r,D)=1$. The other row is
$$
                      C=(ue r,ud r,(2de-u(d+e))r+1).
$$
All three entries are positive and divisible by $u$.

At each weight prime the full power divides either all of $A$ or all of $C$, since $\gcd(u,V)=1$. Every odd predecessor power belongs wholly to $d$, $e$ or $r$. In $C$ the residues are $1,0,0$ up to permutation at all those powers. For example modulo $d$ the first residue is one and the second zero, fixing the third by the sum $N\equiv1\pmod d$. Modulo $r$, only the third part has residue one. The row $A$ avoids 2 at every layer. These are all the required lowest-layer assignments.

## Prescribing the exact valuation

For an odd prime $p\le B$ omitted from $m$, require $uV\equiv-1\pmod p$. This excludes $p$ from the endpoint coefficient. For $p^j\parallel m$, choose $b$ minimal with $p^b>u$ and require
$$
                     V\equiv-p^j\pmod{p^{b+2j}} .
$$
Together with the conditions at $u,D^2$, these are compatible by the Chinese remainder theorem and specify one residue modulo an odd integer $L$.

Now $v_p(N)=j$ and $p\nmid N-1$. The first part of $C$ is a $p$-unit, as $p\nmid ude$ and $2de r=N-1$. The factorial decrement identity implies $p^j\mid K_C$. In $A$ drop $j$ initial zero places. At each layer $w=p^h$ with $b\le h<b+j$ the three residues after this shift are $w-1,w-z,w-z^2$. Their sum $3w-u$ exceeds $w$. Each layer contributes at least one to $v_p(K_A)$ by Legendre's formula. Thus $p^j$ divides all three coefficients, and the endpoint has valuation exactly $j$. Binary avoidance excludes 2. For $m=1$, apply the construction prescribing one additional prime $q>B$ to the first power, with cutoff enlarged to $q$. This guarantees $G>1$.

### Including a prescribed power of two

If $2^h\parallel m$ with $h\ge1$, use the same construction with these changes. Choose $a$ a power of three at least $\max(B,h+1)$ and less than three times that quantity. The proof that $u$ has no prime divisor through $B$ still applies. Replace $2D$ in $N-1=2Dr$ and in the congruence for $uV$ modulo $D^2$ by $2^hD$. Use
$$
                     C=(ue r,ud r,(2^hde-u(d+e))r+1).
$$
Keep the odd-prime conditions for the odd part of $m$. Start $V$ with the fixed binary digits of $1+2^h$ and no other ones through position $h$. Then $u\equiv1\pmod{2^{h+1}}$, $v_2(N-1)=h$ and $r$ is odd. All three parts of $C$ are odd, so their residues already sum to more than one at every layer through $2^h$. This gives $v_2(K_C)\ge h$.

Before starting the free binary word, put $h$ separate pairs of ones into $V$ at positions $t_i,t_i+a$. Use positions above $h+2a$ with the starting positions more than $3a$ apart. At each position $t_i+a$, both $V$ and $zV$ have a one. Each such overlap forces a carry to the next place, giving $v_2(K_A)\ge h$. The lowest-layer assignment at $2^h$ still holds on $A$, whose residues are $1,0,0$. The endpoint valuation is exactly $h$. All these choices consume $O_m(a)$ positions. Leave a gap above the fixed prefix and use the same free word and residue correction as below. This proves the extension to even $m$, including powers of two.

## Counting the weights

In a free $k$-bit word, positions in each residue class modulo $a$ constitute a chain. Within each chain, selected positions must be separated by at least three steps. The number $F_j$ of allowed words on a chain of length $j$ has $F_0=1,F_1=2,F_2=3$ and $F_j=F_{j-1}+F_{j-3}$ for $j\ge3$, by splitting off the last bit. Hence $F_j\ge\lambda^j$. The $a$ chains are independent, giving at least $\lambda^k$ free words.

Any free word extends injectively to an odd admissible $V$ in the required residue class modulo $L$ using at most $k+O(aL+a)$ bits. Put a one at position zero and place the free word starting at $2a+1$. Above it leave $2a$ zero positions and allow $L-1$ further positions separated by $2a+1$. Each available bit gives a unit in the additive group modulo $L$, so its translates generate that whole group. Starting with the singleton residue zero, adjoining the option to use a bit increases the covered set by at least one unless it already fills the group. Otherwise the previous nonempty set would be invariant under translation by a generator. Thus the high positions can repair any residue, and all the required bit separations are preserved. Taking $k=\log_2 X-O_{m,B}(1)$ proves the fixed-cutoff bound.

For the growing cutoff, $\log u=O(B)$ and $\prod_{p\le B}p\le\exp(O(B))$. For the latter bound, the product of primes in $(n,2n]$ divides $\binom{2n}{n}\le4^n$, and one can sum over dyadic intervals. Thus $\log d,\log e=O(B)$ for the displayed choices. Each extra prescribed modulus has logarithm at most $\log u+(2j+1)\log p$. Consequently $\log L=O_m(B)$ for fixed $m>1$, including the modification at 2. For $m=1$ choose the auxiliary prime in $(B,2B]$ by Bertrand's postulate. The same bound holds with an absolute constant since just one valuation is prescribed. A sufficiently small $c_m$ makes the bit overhead, including $\lceil\log_2 u\rceil$ for the weight, $O_m(\sqrt{\log X})$. The claimed counting bound follows.

## An explicit binary test and its limits

Take any odd $V\equiv234928\pmod{455175}$ whose binary ones have gaps of at least three positions. Form $N=7V$, $r=(N-1)/510$ and
$$
              A=(V,2V,4V),\qquad C=(238r,210r,62r+1).
$$
Here $455175=7\cdot255^2$. The congruence and oddness give $V\equiv39853\pmod{65025}$, $V\equiv1\pmod7$, $r\equiv37\pmod{255}$ and $r\equiv1\pmod7$. The binary condition gives $V\equiv1\pmod8$, so $2\parallel N-1$. All the odd full powers in $N-1$ belong to the three pairwise coprime groups $15,17,r$. The first part of $C$ is one modulo 15 and the second is one modulo 17. Both are multiples of $r$. Their complementary zero residues follow from $15\mid210$ and $17\mid238$. At the weight powers $A$ takes $V$ and $C$ takes 7. After dividing $C$ by 7, the first two digits modulo 7 are 6 and 2, forcing a carry. Thus 7 divides both coefficients and the endpoint while $A$ is odd.

There is such a pair at every **$N=8^{99+4760k}-1$, $k\ge0$**. The factor $V=N/7$ has ones exactly three positions apart, and $V\equiv99\equiv1\pmod7$ by the geometric sum. The congruences $4\cdot8^{98}\equiv9436$ and $8^{4760}\equiv1\pmod{65025}$ give $r\equiv37\pmod{255}$. [Lean checks the modular progression and the base-seven obstruction](BinaryLayerWitness.lean).

More generally, **each interval $[7\cdot2^h,7\cdot2^{h+1})$ for $h\ge48$ contains at least $\lambda^{h-48}$ such weights**, where $\lambda=\lambda_3$. A 41-bit word with ones at least three positions apart realizes every residue modulo $455175$. To verify this, start $R_{-2}=R_{-1}=R_0=\{0\}$ and compute $R_j=R_{j-1}\cup(2^{j-1}+R_{j-3})$ in the residue group. At $j=40,41$ the sizes are 455,152 and 455,175. For $k=h-48$, put
$$
                     V=1+8w+2^{k+5}z+2^{k+48}.
$$
Choose any allowed $k$-bit word $w$ and use the 41-bit $z$ to repair the residue, possible since $2$ is a unit. Adjacent blocks keep the required gaps. Here $F_k\ge\lambda^k$ from $F_0=1,F_1=2,F_2=3$.

For a smaller example at $N=245623$, take $(35089,70178,140356)$ and $(25662,109746,110215)$. The former is binary carry-free, and the full lowest-layer assignments all pass. Direct factorial arithmetic gives common endpoint primes exactly $3,7,13,47,67$. This example and the families above isolate the need for higher carries at odd primes in any proof using the [two-row reduction](all_weights_endpoint_reduction.md).

Run `python3 check_prescribed_common_prime.py` and `lean BinaryLayerWitness.lean`. The Python checks use exact coefficients on the small example, separated binary representations and modular inverses for the general construction, and two independent computations of the 41-bit covering.
