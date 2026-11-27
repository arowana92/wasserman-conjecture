# Square and cube exponents for every odd cofactor through 65,535

Put $N=1+d2^a$, with $d$ odd and $a\ge1$. Wasserman's full three-row conjecture holds in each of the following cases.

- $d\le65535$ and $a$ is a square or a cube. More generally $a=x^m$ for any positive $x$ and any degree $m$ divisible by one of $2,3,245,533,637,1435$.
- $d=5Q$, where $Q<2^{23}$ is any odd prime power, and $a=x^m$ with $m$ divisible by one of $2,3,5,11,1283$.
- At every exponent for every odd $d\le65535$ except possibly $12031,56279,60829,61423,65095,65243$. Thus 32,762 of the 32,768 cofactors work without an exponent restriction. Among proper perfect powers $a=x^m$, $m\ge2$, only $12031,56279,65243$ remain. When $a=\ell^m$ for a prime $\ell$, only $12031$ remains and $\gcd(m,2\cdot3\cdot7\cdot13)=1$ is necessary.

The starting point is the exhaustive binary cofactor classification and its refined exponent congruences. A putative counterexample has a row $A=(xq,yq,zq+1)$ with $q=2^a$ and $x+y+z=d$. The classification gives eleven rows and specifies the primes that each coefficient of $A$ must avoid. Three of the exclusions below also follow by applying the earlier independent result for products of two odd prime powers through 1,000. [Cofactor classification](binary_cofactor_rays.md) · [Refined exponent classes](binary_exponent_progressions.md) · [Two-factor calculation](binary_two_factors.md)

## Four exclusions at all exponents

Write $[b]_w$ for the least nonnegative residue of $b$ modulo $w$. If the row avoids $p$, then its three residues sum to $[N]_{p^j}$ at every layer. In particular no individual residue exceeds $[N]_{p^j}$.

At each row in the table, both $p\mid d$ and $p\mid s$, where $sq$ is the indicated unshifted part. Knowing $q$ modulo $w/p$ therefore fixes both entries in the last two columns, for every lift.

| $d$ | $p$ | $s$ | $w$ | $q\bmod(w/p)$ | $[sq]_w$ | $[N]_w$ |
|---:|---:|---:|---:|---:|---:|---:|
| 33661 | 821 | 32840 | $821^2$ | 820 | 641201 | 640381 |
| 62463 | 443 | 53603 | $443^2$ | 440 | 35440 | 8861 |
| 62809 | 587 | 42264 | $587^2$ | 586 | 302305 | 281761 |
| 15871 | 59 | 2065 | $59^3$ | 88 | 181720 | 164375 |
| 15871 | 59 | 2065 | $59^3$ | 1858 | 139948 | 119122 |

In the first three rows the exponent congruence fixes the displayed $q\bmod p$. In the last two, $q\equiv29\pmod{59}$. Testing the 59 lifts against the three slopes $(2065,13450,356)$ at $59^2$ leaves exactly 88 and 1858. In each line the part exceeds the total modulo $w$, so $p$ divides the coefficient of $A$. This contradicts its required avoidance. In particular the only previously remaining cofactor for square and cube exponents, 62809, is excluded for every exponent.

## Degrees of perfect powers

At 61423, every candidate $a$ is 2 modulo 4. At 65095 it is 15 modulo 25. Neither can be a proper perfect power. At 60829 a square is excluded modulo 29. The refined class $a\equiv35049\pmod{627270}$ gives $7\mid a$. At degree $m\ge3$ this would give $343\mid a$. The exponent class with this added condition forces a carry in $A$ through $7^4$.

The remaining three rows give these sufficient degree divisors for exclusion. At 56279 the conditions work on both candidate rows.

| $d$ | Degree divisors excluded | Corresponding modulus for the exponent |
|---:|:---|:---|
| 12031 | $2,3,7,13$ | $4,13,113,53$ |
| 56279 | $2,3,41,49$ | $7,7,83,7^4$ |
| 65243 | $2,3,5,13$ | $5,13,41,53^2$ |

At the single-prime entries and at modulus 4, the exponent class contains no residue of the listed power. At 12031 use $a\equiv769251571\pmod{2043279108}$. For the entries $7^4$ and $53^2$, test the candidate row at one additional prime-power layer as follows. Given a residue $b=x^m\bmod p^j$, use CRT with the row's progression $a\equiv a_0\pmod H$. On each compatible class the scalar $2^a\bmod p^{j+1}$ is fixed since $2^{\operatorname{lcm}(H,p^j)}\equiv1\pmod{p^{j+1}}$. In every case the sum of the three part residues overflows at some layer through $p^{j+1}$. The same finite test with $b=0$ gives the assertion at 60829. Taking least common multiples of one excluded degree from each row gives the six sufficient degrees in the first assertion after removing multiples of smaller ones. The previous prime-power exponent classification already excluded every cofactor except 12031 for prime bases $\ell$. [Prime-power exponents](binary_cofactor_rays.md#exponents-that-are-multiples-or-powers)

For $d=5Q$ the earlier calculation leaves two classes of exponents. One is 2 modulo 4, which excludes every proper perfect power. The other has $a\equiv2337019\pmod{2935556}$ and slopes $(4403337,2918400,17158)$. Squares are excluded modulo 4, cubes modulo 13, and eleventh and 1,283rd powers modulo 56,453. For fifth powers test the 125 possible bases modulo 125 by the preceding CRT calculation at $5^4$. Every compatible exponent gives a quinary carry. [Two remaining factor-five rows](factor_five_exponents.md)

Run `python3 check_binary_perfect_powers.py` from the repository root. It exhausts all scalar lifts used in the four unconditional exclusions and checks the modular power tests with integer factorial floors.
