# Proper prime powers at bounded cofactors of $10P+7$

Let $P=p^k$, $R=q^e$ be prime powers of positive exponent with $10P+7=DR$, $1\le D\le10{,}000$. At the weight $N=(2P+1)(5P+1)$ the following results concern the full three-row conjecture.

* **Only finitely many possible exceptions have $k>1$ or $e>1$**, even with both prime bases unrestricted. The proof of finiteness does not give an explicit largest exception.
* **The conjecture holds at every exponent if $p\le10^6$ or $q\le10^6$.** This includes the entire range of $D$.
* For arbitrary prime bases, **$k\ge 300{,}000{,}000{,}000{,}000$ or $e$ at least this large suffices**. The remaining pairs of ordinary primes are still subject to the [logarithmic counting bound](seven_quadratic_suffixes.md). The [cofactor theorem](seven_quadratic_cofactors.md) already settles every input outside three cofactors.

## Exact limits of the binary words

Keep the notation $j,L,n$ from the binary count, and let $h$ be the bit length of $n$. Regard its leading-to-trailing digit graph as a deterministic graph accepting every word that reaches a nonempty subset. Repeatedly delete nodes with no successor in the current graph. Paths in the remaining graph extend indefinitely. At $D=4587,9453$ no pruning is required. At $D=7383$ the process stabilizes after six deletions of layers. Thus removing at most the last six digits from an accepted word leaves an extendible prefix. For a word of length at most six the bound below holds directly.

The verifier follows every live branch from the leading 1 until a node repeats. It checks that each node on the cycle has exactly one live successor. Consequently every infinite word so obtained is eventually periodic. A prefix of length $h_0$ with value $b$ followed by a repeating block of length $v$ and value $c$ has value $(b+c/(2^v-1))/2^{h_0}$ as a binary fraction. This gives exactly the following fractions $f$.

| $D$ | Binary limits $f$ |
|---:|---|
| 4,587 | $256/417$ |
| 7,383 | $19597/35328$, $486233/876544$, $836139/1507328$, $46678369/84148224$ |
| 9,453 | $256/411$ |

At every accepted word $|n-f2^h|\le64$ for some listed $f$, by agreement through all but at most six trailing digits. Since $n=\lfloor(20P^2+24P+7)/(32D2^j)\rfloor$ and $2^j\le P$, this implies
$$
     |P^2-C2^{j+h}|\le1{,}100{,}000P,\qquad C=8Df/5 .
$$
Indeed the numerator error after multiplication by 20 is at most $32D\cdot65P+24P+7$. For $R=(10P+7)/D$, expansion of its square gives the same type of inequality with $C'=100C/D^2$. In particular, writing $Q$ for either $P$ or $R$, at any possible failure there is an integer $m=j+h$ and one of the corresponding constants $c$ such that
$$
                          |Q^2-c2^m|\le 2\cdot10^6 Q.                 \tag{1}
$$
For the $R$ bound use $P<DR/10$ and $D\ge4587$. Each of the twelve positive constants $c$ is between $1/1000$ and $10^4$, its reduced numerator and denominator are below $10^{11}$, and at least two different odd primes occur in its numerator and denominator combined. These are integer checks on the displayed fractions.

## One exponent bound independent of both bases

Suppose $Q=s^t$ with $s$ an odd prime. The [finite theorem](seven_quadratic_suffixes.md) lets us assume $P>2^{65536}$, so $Q>2^{65526}$. The real linear form
$$
                     \Lambda=m\log2+\log c-2t\log s
$$
is nonzero: the odd prime support of $c$ contains at least two primes. Inequality (1) and $|\log(1+x)|\le2|x|$ for $|x|<1/2$ yield $|\Lambda|\le 4\cdot10^6/Q$.

Apply Matveev's Corollary 2.3 to $2,c,s$ over the rational field, with modified heights $A_1=\log2$, $A_2=26$, $A_3=\log s$. Its weighted coefficient parameter is the maximum of $1,mA_1/A_3,A_2/A_3,2t$. Because $|m\log2-2t\log s|\le|\log c|+1<11$, this maximum is at most $30t$. The constant $C_1(3)$ in (2.6) is at most $2^{38}$. It follows that
$$
          t < 16+26\cdot2^{38}(1+\log(30t)).
$$
At $t=3\cdot10^{14}$ the right side is less than $16+26\cdot2^{38}\cdot38<3\cdot10^{14}$. Its derivative is less than one from this point onwards. This proves the uniform exponent bound. The weighted parameter is essential here because $m$ grows with $\log s$. [Matveev, Corollary 2.3 and equation (1.3)](https://doi.org/10.1070/IM2000v064n06ABEH000314)

## Finiteness for higher powers

The three remaining cofactors force both exponents odd. Therefore it suffices to consider $3\le t<3\cdot10^{14}$. Fix $c,t$, and fix the residue $r\in\{0,\ldots,2t-1\}$ of $m$ modulo $2t$. Write $m=2t u+r$ and $\alpha=(c2^r)^{1/(2t)}$. Taking the difference of squares in (1) gives $|s^t-(\alpha2^u)^t|\le2\cdot10^6$. Factoring the difference of $t$-th powers then gives
$$
                \left|s/2^u-\alpha\right|\ll_{c,t,r} (2^u)^{-t}.
$$
For irrational $\alpha$, Roth's theorem applied, for example, with approximation exponent $5/2$ permits only finitely many reduced fractions of this kind. For rational $\alpha=A/B$ in lowest terms, a distinct fraction differs by at least $1/(B2^u)$, and equality fixes its reduced denominator since $s$ is odd. Bounded $u$ also gives bounded $s$. There are finitely many $c,t,r$, so all sufficiently large possible exceptions in the entire cofactor range have both $P,R$ ordinary primes. [Roth, *Rational approximations to algebraic numbers*](https://doi.org/10.1112/S0025579300000644)

## All exponents at prime bases through one million

Put $\theta=2\log_2 s$, $\mu=\log_2 c$. The upper bound on $\Lambda$ requires $|t\theta-m-\mu|<10^7/Q$. Since $Q>2^{65526}$, it suffices to exclude even a much larger tolerance. Write $\|x\|$ for distance to the nearest integer. For each base and offset it suffices to exhibit an integer $v<10^{32}$ with
$$
                   \|v\mu\|-10^{18}\|v\theta\|>1/1000 .
$$
The triangle inequality then gives $|t\theta-m-\mu|>1/(1000v)$ for every $t<10^{18}$, including every exponent allowed by Matveev.

The first program uses directed integer squaring for binary logarithms to 200 bits and tries convergent denominators of the resulting rational approximation to $\theta$. Given $A=\lfloor2^{200}\theta\rfloor$ and $B=\lfloor2^{200}\mu\rfloor$, the uncertainties in each of $v\theta$ and $v\mu$ are at most $v/2^{200}$. All success tests include these errors. The choice of convergents speeds the search; the final certificate is the displayed integer-distance inequality.

The second program uses $\theta'=\log_2 s$ and coefficient $2t<6\cdot10^{14}$. It encloses natural logarithms by the positive series for $2\operatorname{arctanh}((x-1)/(x+1))$, with integer rounding in both directions and a geometric tail bound. It then searches denominators independently, requiring a distance surplus greater than $1/256$. Both cover all 78,497 odd prime bases and all twelve offsets, or 941,964 checks. The largest final denominators used were 60,094,944,827,760,486,551,335 and 2,440,480,511,603,505,239,330, respectively. The binary base is already excluded for both powers in any surviving shape.

```sh
python3 check_seven_quadratic_bases.py 1000000
python3 check_seven_quadratic_bases_cross.py 1000000
```
