# The exact pair-gcd bound on a homogeneous plane

Fix positive $r,s,t$ with $\gcd(r,s,t)=1$. On the plane $rx+sy=tz$ write $u=x+y+z$ and
$$
 H=\operatorname{lcm}(\gcd(y,z),\gcd(x,z),\gcd(x,y)).
$$
The ratio $u/H$ is invariant under scaling. Its minimum $\mu(r,s,t)$ is attained, and the following formula computes it by a finite search. For each $d\mid r$, $e\mid s$ with $\gcd(d,e)=1$, let $f$ be the largest divisor of $t$ coprime to $de$, and $q=t/f$. Minimize
$$
       \left(1+\frac rt\right)\frac{i}{d}
       +\left(1+\frac st\right)\frac{k}{e}                            \tag{1}
$$
over these choices and all $1\le i,k\le q$ with $(r/d)i+(s/e)k\equiv0\bmod q$.

It follows from the [companion divisibility theorem](ratio_companion_divisor.md) that a row on this plane and any opposite row with part ratio $a:b$ share an endpoint prime at every weight whenever $C(a,b)/(a+b)\le\mu(r,s,t)$. The minimum is exact for the direction inequality used here.

To prove the formula, start with a primitive row. Its three pair gcds $d,e,f_0$ are pairwise coprime divisors of $r,s,t$, respectively. Write its entries as $(i e f_0,k d f_0,l d e)$. Then $(r/d)i+(s/e)k=(t/f_0)l$. Its value of $u/H$ equals (1). Since $f_0\mid f$, replacing the modulus $t/f_0$ by $q$ preserves the congruence. Replacing each of $i,k$ by its least positive residue modulo $q$ can only decrease (1).

Conversely take any admissible choices in the finite search. With $l=((r/d)i+(s/e)k)/q$, the positive row $(i e f,k d f,l d e)$ is on the plane. Its pair-gcd lcm is a multiple of $def$, so its actual $u/H$ is at most (1). Taking the minimum on both sides proves equality, including attainment by a primitive row after normalization.

## Coefficients supported on few primes

If the three coefficients admit **no three distinct primes, one dividing each coefficient**, then $\mu>1$. Otherwise a row with $H\ge u$ would have three pair gcds greater than one. Indeed, when one pair gcd equals one, the other two have product dividing the coordinate at that position, so $H<u$. Choosing a prime from each of three pairwise coprime nontrivial pair gcds would give the forbidden distinct primes.

Thus any such plane works against every integer ratio and every [geometric-sum ratio](part_ratio_predecessor.md), at every weight. Examples include any plane with a unit coefficient, any plane whose coefficients together use at most two primes, and any plane with two coefficients that are powers of the same prime. Also, **two planes both having this property work against each other**, even with unrelated coefficient primes. A counterexample would require coprime primitive sums $u,v$ dividing $N$ and $N-1\mid\operatorname{lcm}(H_1,H_2)$, whereas $H_1H_2\le(u-1)(v-1)<N-1$.

More generally, if two planes have $\mu_1\mu_2>1$, a failure requires $N\le\mu_1\mu_2/(\mu_1\mu_2-1)$. This follows from $N\ge uv\ge\mu_1\mu_2 H_1H_2\ge\mu_1\mu_2(N-1)$.

## A closed form

If $t$ is coprime to $rs$ then $q=1$ throughout and $i=k=1$ suffices. In particular for every $\gcd(R,t)=1$,
$$
 \mu(R,R,t)=\left(1+\frac Rt\right)
       \min_{\substack{de=R\\ \gcd(d,e)=1}}\left(\frac1d+\frac1e\right). \tag{2}
$$
Any coprime pair of divisors of $R$ can be enlarged to a coprime factorization of $R$, decreasing or preserving the cost. As an example $\mu(6,6,31)=185/186$, attained at $(93,62,30)$. This proves the endpoint prime against every integer ratio $1:b$ with $b\le185$ on this plane.

Among the 27,187 primitive normals $1\le r\le s\le40$, $1\le t\le40$, exactly 7,887 have $\mu>1$. Of these, 4,867 already have no distinct prime representatives. The [verifier](check_plane_gcd_ratio.py) checks each attained ratio against its primitive witness. An independent enumeration of every possible direction with $H>u$ checks the threshold and all minima below one for every one of the 171 normals through height seven. Equality $H=u$ is impossible for primitive triples since $\gcd(H,u)=1$.

~~~sh
python3 check_plane_gcd_ratio.py
python3 check_plane_gcd_ratio.py --plane 6 6 31
~~~
