# A cube of 1,024 exponents at each of the bases 2, 3 and 5

**Theorem.** The full three-row conjecture holds at every
$$
             1+5^a 2^b 3^c,\qquad 1\le a,b,c\le1024 .
$$
Two independent integer programs verify all 1,073,741,824 exponent triples. For each weight, the row forced by a counterexample would have to avoid $2,3,5$. A single binomial factor of its coefficient is divisible by at least one of them.

**Translated cubes.** Let $A$ be any nonnegative multiple of $2^{1086}3^{1055}$. The whole cube also works with $(a,b,c)$ replaced by $(a+A,b+B,c+C)$ whenever the exponents stay positive and
$$
   2\cdot3^{31}\mid B,\qquad 2^{62}\mid C,\qquad
            2^B3^C\equiv1\pmod {5^{A+1024}}. \tag{1}
$$
Negative shifts in (1) use modular inverses. For each $A$, the shifts $B,C$ can independently be arbitrary nonnegative multiples of $4\cdot3^{31}5^{A+1023}$ and $2^{62}5^{A+1023}$, respectively.

At $A=0$, the two-dimensional lattice (1) has a basis of positive vectors whose coordinates have at most 375 decimal digits. **Every square of $10^{375}$ consecutive positive exponents on each of the $b,c$ axes contains a translated $1024$-by-$1024$ square that works simultaneously for all $1\le a\le1024$.**

## Checking the original cube

Put $D=5^a$, $Q=2^b$, $R=3^c$, $W=QR$. In a counterexample the [forced row](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) avoids $2,3,5$ and contains the part $xW$ with $x=W^{-1}\bmod D$ in $[1,D-1]$. The [endpoint theorem](all_weights_endpoint_reduction.md) applies to arbitrary triples of rows. By grouping the other two parts, this row must pass all three binomial carry tests. Write $u\preceq_s v$ when every base-$s$ digit of $u$ is at most the corresponding digit of $v$. The tests inspect the specified low digits of each side.

| Base | Required digit containment | Digits inspected | Triples left |
| ---: | --- | ---: | ---: |
| 2 | $xR\preceq_2 DR$ | 64 | 3,182 |
| 3 | $xQ\preceq_3 DQ$ | 32 | 7 |
| 5 | $(xW-1)/D\preceq_5 W$ | at most 3 on the remaining triples | 0 |

At base 2 the low $b$ digits of the cut contain just the units digit in its complement. The high digits give the binary test. The base-3 test follows by removing $c$ digits. At base 5 the low $a$ places of both $xW$ and $1+DW$ equal one, giving the last comparison. This uses the multinomial carry criterion. [Bergman, Lemma 5](https://arxiv.org/html/0806.0607)

Here $j$ is the first layer at which the residue of $(xW-1)/D$ modulo $5^j$ exceeds the residue of $W$. Every survivor of the first two tests is listed.

| $a$ | $b$ | $c$ | $j$ |
| ---: | ---: | ---: | ---: |
| 1 | 4 | 1 | 1 |
| 1 | 24 | 1 | 1 |
| 1 | 4 | 2 | 3 |
| 1 | 6 | 2 | 2 |
| 2 | 1 | 3 | 2 |
| 2 | 7 | 5 | 2 |
| 6 | 4 | 4 | 3 |

The Python verifier fixes $a,c$ and walks backward in $b$. Decreasing $b$ by one doubles $xR$ modulo $DR$. It tests the low binary bits, then factorial-floor differences in bases 3 and 5. The independent C++ verifier fixes $b,c$ and reads the digits of $W^{-1}$ in base 5 to recover every $x$ as $a$ increases. If the next 24 base-5 digits of the inverse above $x$ form $z$, the quotient $(xW-1)/D$ is congruent to $-zW\bmod 5^{24}$. This supplies the last test without forming the quotient. [Python floor search](check_exponent_congruence_cube.py) · [Independent base-5 search](check_exponent_congruence_cube_cross.cpp)

## Moving the 5-exponent

Set $D'=5^{A+a}$. The choice of $A$ gives
$$
                 D'\equiv D\pmod{2^{1088}3^{1056}}.
$$
For binary powers this follows by repeated squaring from $5\equiv1\bmod4$. At the ternary power it follows from Euler's theorem. Let $t=(xW-1)/D$, so $1\le t<W$. The inverse of the **original** $W$ modulo $D'$ is
$$
              x'=\frac{1+tD'}{W}
                  = x+t\frac{D'-D}{W}.
$$
Indeed $W$ divides $D'-D$, the displayed value is positive and below $D'$, and $Wx'=1+tD'$. Since $W=2^b3^c$ with $b,c\le1024$, we obtain $x'\equiv x\bmod 2^{64}3^{32}$. The base-5 quotient $t$ is exactly the same after increasing $a$. This step is what allows the $a$-exponents of the translated cubes to be arbitrarily large.

Now replace $Q,R,W$ by $Q',R',W'$ after shifting $b,c$. Condition (1) implies $Q'\equiv Q\bmod3^{32}$, $R'\equiv R\bmod2^{64}$ and $W'\equiv W\bmod5^{A+1024}$. The new inverse modulo $D'$ is still $x'$. Thus both grouped tests have their original digits at 2 and 3. For the seven triples that reach base 5, $a+j\le9\le1024$. Consequently $(x'W'-1)/D'\equiv t\bmod5^j$ and $W'\equiv W\bmod5^j$. The witness at 5 transfers as well. Lean 4.36.0 checks the two integral identities for the shifted inverse and quotient, and the bound placing the new inverse in its required range. [Formal arithmetic lemmas](ExponentCubeTransfer.lean)

## Periods in the two other exponents

For $A=0$, put $m=5^{1023}$. The order of 2 modulo $5m$ is $4m$, since 2 generates the units modulo 5 and $v_5(2^4-1)=1$. Let $\ell$ be the discrete logarithm of 3 to base 2 modulo $5m$. The binary shift condition and reduction of the last condition of (1) modulo 5 force $4\mid B$. Write $B=Uu$, $C=Vv$ with $U=4\cdot3^{31}$ and $V=2^{62}$. The last condition becomes
$$
                 3^{31}u+2^{60}\ell v\equiv0\pmod m.
$$
It describes an index-$UVm$ lattice in the original shift coordinates. The verifier lifts the logarithm one base-5 digit at a time, constructs a basis, and reduces it by the Euclidean algorithm for a two-dimensional lattice. It outputs two positive vectors $v_1,v_2$, checks their congruences directly by modular powering, and checks $|\det(v_1,v_2)|=UVm$. It also verifies
$$
       1024+\max_{j=1,2}(|(v_1)_j|+|(v_2)_j|)<10^{375}.
$$
Express the difference between the center of the desired square and the center of $[1,1024]^2$ in this basis and round both coefficients to nearest integers. In each coordinate the rounding error is at most half the corresponding sum of absolute coordinates. This places the translated block inside the desired square. [Lattice construction and arithmetic checks](check_exponent_orbits.py)

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 -Wall -Wextra check_exponent_congruence_cube_cross.cpp -o .work/exponent_cube
.work/exponent_cube
python3 check_exponent_congruence_cube.py
python3 check_exponent_orbits.py
lean +leanprover/lean4:v4.36.0 ExponentCubeTransfer.lean
~~~
