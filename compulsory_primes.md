# Primes that every certificate must contain

A certificate at $N$ is a set of primes containing a common divisor for every choice of three positive trinomial rows. A selection whose coefficients have gcd exactly $p$ makes $p$ compulsory in every certificate, regardless of its size.

* At 785 there are **33 minimum certificates**, each of size four. Their intersection is empty.
* At 699 the compulsory primes are exactly **233 and 349**. At 717 they are exactly **179 and 239**. Both weights have minimum certificate size five.
* For every $0<\epsilon<1/4$, at $\gg_\epsilon X/(\log X)^2$ composite weights in $[X,2X]$, the full conjecture holds, no prime divides every coefficient, and every certificate contains a prime at least $X^{1/4-\epsilon}$.
* At an explicit 31-digit weight $N$, **both** $N/3$ and $(N-1)/2$ are compulsory primes. There are minimum certificates of size four using either 2 or 3.

## Prescribing the exact gcd

Write $K_R$ for the coefficient of row $R$. Fix coprime integers $d\ge1$, $h\ge2$. There exist fixed positive $a,b$ and an effectively constructible reduced arithmetic progression of primes $q$ such that for every sufficiently large prime in that progression, $N=hq+1$ is a multiple of $d$ and
$$
\gcd(K_E,K_A,K_C)=N/d,\qquad
E=(1,1,N-2),\quad A=(a,b,N-a-b),\quad C=(1,q,(h-1)q).
$$

For each $p^e\parallel d$ impose $a\equiv d\bmod p^{e+1}$, and for each $p^e\parallel h$ impose $a\equiv1+p^e\bmod p^{e+1}$. Choose a positive CRT solution. Let $b$ be a product of powers of the primes dividing $dh$, with the power of each prime larger than $\max(a,dh)$ and divisible by $p^{v_p(dh)+1}$. Thus $\gcd(a,b)=d$, and $a,b$ occupy disjoint ranges of digit positions in each of these bases. For each $p\mid dh$ take $p^{j_p}>a+b+dh$ and prescribe
$$
hq+1\equiv a+b\pmod {p^{j_p}}.
$$
If $p^e\parallel h$, then $p^e\parallel a+b-1$. After cancelling $p^e$ the congruence fixes a unit class for $q$ modulo $p^{j_p-e}$. If $p\mid d$, its right side after subtracting one is a unit and $h$ is invertible. The resulting progression for $q$ is reduced. It contains infinitely many primes by Dirichlet's theorem.

For $N>a+b$, the third entry of $A$ is a multiple of $p^{j_p}$ at every prescribed base, so $A$ avoids all primes of $dh$. Also $d\mid N$ and $\gcd(d,N/d)=1$. The integer $K_A A_i/N$ is the coefficient with the $i$-th part and the weight each reduced by one. Hence $N/\gcd(N,a,b)=N/d$ divides $K_A$. Both $K_E$ and $K_C$ are multiples of $N$ because their rows contain 1. The endpoint restricts the gcd to divisors of $N(N-1)=d(N/d)hq$, and $C$ avoids $q$ for $q>h$ by its two base-$q$ digits. These facts give the asserted equality.

If $d\ge3$ and $P=N/d$ is also prime with $P>d$, then for $q>\max(a,b,h)$ a second exact gcd is
$$
\gcd(K_E,K_A,K_D)=q,\qquad D=(P,P,(d-2)P).
$$
Indeed $D$ avoids $P$ by its one occupied base-$P$ column. Each of $A,D$ has two entries nonzero modulo $q$, which forces a carry at weight $1\bmod q$. The endpoint has just one factor $q$. Together the two gcds make $P,q$ compulsory whenever the weight has a certificate. The explicit examples below realize both prime gcds.

## A family with a lower bound on the size of the primes

Take $h=2$, $d=3$ and primes $q\equiv7\bmod216$, with $N=2q+1>15$. Then $N\equiv15\bmod432$. One can take $a=3$, $b=12$. Indeed $3,12$ have disjoint binary digits and add freely in base 3 as $10_3+110_3=120_3$. The third part is divisible by both 16 and 27. The preceding proof gives gcd $v=N/3$.

The [endpoint reduction and two-factor predecessor theorem](all_weights_endpoint_reduction.md) prove the full conjecture at every such weight. There is no universal prime. An odd base has odd digit sum at odd $N$, and a digit sum of one would make $N$ a prime power, impossible since $3\parallel N$ and $N>15$. In base two, $N\equiv15\bmod16$ already has four ones. At any prime with at least three digit units these units can be divided into three nonempty parts without a carry.

We sieve the integers $v=(2q+1)/3$ for the primes in this progression with $X\le2q+1\le 2X$. Their prime-counting main term $Y$ is the logarithmic-integral mass for $q$ divided by $\varphi(216)$, so $Y\asymp X/\log X$. Every $v$ is coprime to 6. For squarefree $s$ coprime to 6, $s\mid v$ is the single reduced class $2q\equiv-1\bmod s$, combined with the original class by CRT. Thus the local main term is $Y/\varphi(s)$. For any fixed $\alpha<1/2$ the sum of the absolute errors over $s\le X^\alpha$ is $O_{\alpha,B}(X/(\log X)^B)$ for every fixed $B$. Apply Bombieri–Vinogradov at the two endpoints to the distinct moduli $216s$. [Bombieri, On the large sieve](https://www.cambridge.org/core/journals/mathematika/article/abs/on-the-large-sieve/4AE92826645F6F3393D30AF0D2767346)

The local densities $1/(p-1)$ at $p>3$ and zero at 2 and 3 have sieve dimension one. Indeed the ratio of the local sieve factor to $1-1/p$ is $1-1/(p-1)^2$ at $p>3$, with convergent product. The sieve product through $z$ is therefore of order $1/\log z$ by Mertens' theorem. For $z=X^{1/4-\epsilon}$ choose $1/2-2\epsilon<\alpha<1/2$. The lower linear sieve has strictly positive lower function at $\alpha/(1/4-\epsilon)>2$. It gives $\gg_\epsilon Y/\log z$ values of $v$ with every prime factor at least $z$. Each gives a different certified weight, and the exact gcd $v$ forces every certificate to contain one of those factors. [Matomäki and Zuniga–Alterman, Definition 2.4 and Lemma 2.5](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/weighted-sieves-with-switching/986429394BA969687D48224E177E2F1C)

## Two specific compulsory primes at one weight

Take
$$
P=1000000000000000000000003631621,\qquad
N=3P,\qquad q=(N-1)/2,\qquad r=N-2 .
$$
Here $P,q,r$ are prime and $N\equiv15\bmod432$. Thus the two displayed gcds give $P$ and $q$. Each prime among $P,q,r$ has a unique avoiding row up to permutation, respectively $D=(P,P,P)$, $C=(1,q,q)$ and $E=(1,1,r)$, since its digit units must be split among three nonempty parts. These rows are distinct, and all three coefficients are divisible by both 2 and 3. For 2, the row $D$ repeats an odd part, the row $C$ has two odd $q$'s, and the endpoint is even. For 3, each row has a part coprime to 3 while $3\mid N$.

It follows that both $\{2,P,q,r\}$ and $\{3,P,q,r\}$ certify the weight. In each list no row can avoid two primes. All primes are individually avoidable as above, so the minimum is four. In fact every four-prime certificate contains $P,q$, exactly one of 2 and 3, and one prime outside $N(N-1)$. Without either 2 or 3, $E$ avoids both extra primes, and $C,D$ cover $q,P$. With both 2 and 3, $A$ avoids both and again $C,D$ cover $q,P$. The same conclusions hold at the smaller example $N=306303$.

## The three finite classifications

At 785 a prime in a minimum list must be at most 785. All primes are individually avoidable. A list of four therefore certifies exactly when no row avoids two list members. The endpoint avoids every prime outside $\{2,5,7,157\}$, so at least three list members belong to this set. The table gives all options for the fourth prime from an exhaustive check on all 51,352 unordered positive rows.

| Three specified primes | All possible fourth primes |
|---|---|
| 2, 5, 7 | 3, 19, 23, 29, 41, 71, 157 |
| 2, 5, 157 | 3, 7, 17, 19, 23, 29, 41, 71, 97, 389 |
| 2, 7, 157 | 3, 5, 19, 23, 29, 37, 41, 71, 257, 769, 773 |
| 5, 7, 157 | 2, 3, 11, 19, 23, 29, 41, 71 |

Only the four core primes together occur in more than one row of the table. Deduplication gives 33 lists with empty intersection.

At 699 the rows $(1,1,697),(3,48,648)$ give gcd 233 on adjoining $(1,349,349)$ and gcd 349 on adjoining $(233,233,233)$. At 717 the rows $(1,1,715),(9,192,516)$ give gcd 179 on adjoining $(239,239,239)$ and gcd 239 on adjoining $(1,179,537)$. At each weight the verifier exhibits four five-prime certificates whose intersection consists of the two specified primes. Their [minimality was checked previously](RESEARCH.md). The intersection of all minimum certificates and that of all certificates therefore coincide in each case.

Run the two commands below with python-flint installed for the second one. The verifier uses exact binomial products for the finite classifications and 41 progression weights through 50,000. A separate factorial-floor calculation checks the construction for 78 pairs $(d,h)$, including three with two prime gcds. At the large example it checks recursive integer primality certificates and uses Legendre's factorial floors at every prime of the endpoint coefficient to recover the exact gcds.

~~~sh
python3 check_wasserman_785.py
python3 check_compulsory_primes.py
~~~
