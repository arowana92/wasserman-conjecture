# All outer exponents at every Mersenne prime

Let $R=2^h-1$ be any Mersenne prime. **For every $1\le d\le104$ and every $e\ge1$, Wasserman's three-row conjecture holds at $1+dR2^e$.** For $1\le d\le512$ the same is true for at least 442 multipliers. The test covers all Mersenne primes without a bound on $h$. For example $d=3\cdot5\cdot7$ is the first multiplier the unrestricted-height test leaves open.

Here is the full list of odd multipliers through 511 left by the test. All other odd multipliers in that range work. An even multiplier reduces to its odd part by increasing $e$.

~~~
105 117 161 171 177 187 195 205 209 213 217 219 221 225 231 245 253
273 285 291 309 323 329 345 357 369 371 377 385 403 405 413 425 429
435 437 441 445 451 459 465 469 473 475 477 481 483 485 493 497 505
~~~

The result extends to every $R=2^h-1>1$ that is a prime power. Indeed if $2^h-1=q^\beta$, reduction modulo four shows $\beta$ is odd. For $\beta>1$ the quotient $(q^\beta+1)/(q+1)$ is an odd integer greater than one dividing a power of two. Thus $R$ is prime.

## Four bounded integers

Take $d$ odd and put $H=2^h$, $W=2^e$, $N=1+dRW$. First suppose $H\ge2d$. Then $dR<2^{3h}-1$. By [the adjacent prime-power argument](adjacent_prime_power_weights.md), in a counterexample the nonendpoint rows must assign 2 and $R$ to opposite rows. Their respective weight scales $U,V$ satisfy $UV=N$. Each is at least three and coprime to $dR$. The row avoiding 2 has parts
$$
           WUi,\quad WUj,\quad W(dR-U(i+j))+1.
$$
Here $i,j$ are positive, the third upper part is positive since its actual part is divisible by $U$, and $Ui,Uj$ have disjoint bits contained in $dR$. Order $i<j$. The other row has parts $VRa,VRb,V(U-R(a+b))$ for positive $a\le b$. Hence
$$
                 R(a+b)<U,\qquad U(i+j)<dR,\qquad (i+j)(a+b)<d.
$$
Also $v_2(i)\ne v_2(j)$ since $U$ is odd and the products are disjoint.

At a full prime power $w$ of $d$, one of the rows must have residues $0,0,1$ modulo $w$. Thus $w$ divides one of $i,j,a,b$. If both $i,j$ are zero modulo $w$, or both $a,b$ are, its lowest-layer condition can be dropped. Otherwise assign $w$ to one eligible row. On the first row exactly one of $i,j$ must vanish and the condition is $WU(i+j)\equiv1\bmod w$. On the second row exactly one of $a,b$ must vanish and $U\equiv R(a+b)\bmod w$, since $V\equiv U^{-1}\bmod w$. Discard assignments for which the relevant sum fails to be a unit. Let $A,C$ be the products assigned to the first and second row in this way. Necessary conditions are
$$
 U(i+j)\in\langle2\rangle\bmod A,\qquad U\equiv (H-1)(a+b)\bmod C. \tag{1}
$$
The first condition uses just the existence of the exponent $e$.

## A finite graph for arbitrarily long runs of binary ones

Let $s=a+b$ and $l=\lceil\log_2(2d)\rceil$. For $h\ge l$ write $U=vH+u$, $0\le u<H$. We must have $v\ge s$. Otherwise $U=sH-t$ with $1\le t<s$ by $U>(H-1)s$. Both $it,jt<d\le H/2$. The two residues of $iU,jU\bmod H$ would then both contain the bit $H/2$.

The upper block of $dR=(d-1)H+(H-d)$ is $d-1$. Its lower block has the first $l$ bits of $-d$ followed by arbitrarily many ones. Process the binary digits of $u$ from low to high. At position $k$, put $u_k=u\bmod2^k$ and track
$$
 x=\lfloor iu_k/2^k\rfloor,\quad y=\lfloor ju_k/2^k\rfloor,\quad
 \rho=(u_k+\sigma)2^{-k}\bmod AC,
$$
where $\sigma\equiv0\bmod A$ and $\sigma\equiv s\bmod C$. The initial state is $(0,0,\sigma\bmod AC)$. Read an input bit $b$, taking $b=1$ at position zero because $U$ is odd. The output bits are the parities of $x+ib,y+jb$. Discard the transition if their sum exceeds the target bit in $dR$. The new state is
$$
 \left(\left\lfloor(x+ib)/2\right\rfloor,\,
       \left\lfloor(y+jb)/2\right\rfloor,\,(\rho+b)/2\bmod AC\right).
$$
Modulus one imposes no condition. These transitions follow directly from the three state definitions.

After the first $l$ columns take the complete reachable set under any number of target-one transitions. It is finite because $0\le x<i$, $0\le y<j$, $0\le\rho<AC$. At the boundary a state admits a high part $v\ge s$ only if $iv+x,jv+y$ have disjoint bits contained in $d-1$. Condition (1) requires in addition
$$
       (v+\rho)(i+j)\in\langle2\rangle\bmod A,\qquad v+\rho\equiv s\bmod C.
$$
Indeed $H$ is a power of two, and modulo $C$ we have $v+\rho=(U+s)/H$. All possible $v$ are bounded by $(d-1-x-y)/(i+j)$. Empty reachability excludes the assignment simultaneously for all $h\ge l$ and all $e$. For the optional early high-part test, $x=\lfloor it\rfloor,y=\lfloor jt\rfloor$ for the same $0\le t<1$. At fixed $y$ and $i<j$, only the two endpoint values $\lfloor i y/j\rfloor$ and $\lceil i(y+1)/j\rceil-1$ need be considered.

For $h<l$ and $d\le512$ the only Mersenne primes are $3,7,31,127$, as direct inspection through $h=9$ shows. All odd cofactors $dR$ in these cases are below 65,536. None equals one of the ten cofactors left in the [verified all-exponent classification](binary_cofactor_rays.md), by division using its table. Thus these heights work as well.

## Computation

The standard-library program enumerates 8,299 quadruples and 8,466 assignments for odd $d\le511$. Only 453 assignments at the 51 listed odd multipliers reach an accepting state. For $d\le127$ it also checks the fixed-height version at three successive heights by direct enumeration of odd $U$ and comparison of the actual products. All 381 comparisons agree.

~~~sh
python3 check_mersenne_multipliers.py
~~~
