# Exact divisibility near rational trinomial proportions

Let $s=(s_1,s_2,s_3)$ be a primitive triple of positive integers, $D=\sum s_i$, and $K(x)=(\sum x_i)!/\prod x_i!$. Define $G(s)$ as the largest integer that divides **every** positive trinomial coefficient in some open neighborhood of the proportions $s/D$. The neighborhood is measured in the triangle of proportions. It includes every weight and allows independent changes in the three part sizes.

**For every rational direction, $G(s)$ is exactly computable.** For a prime $p$ and a proposed valuation bound $B\ge0$, a graph on at most $12D(B+1)$ states decides between the following two outcomes.

* Throughout an explicit neighborhood of $s/D$, all coefficients are divisible by $p^{B+1}$.
* There are a fixed nonnegative integer triple $a$ and positive integers $L,h$ such that the rows
  $$
     a+Ls\,(p^{hn}-1)/(p^h-1),\qquad n\ge1,
  $$
  all have valuation at most $B$. The algorithm supplies their shared prefix and repeating digit block.

In particular, near $(1,3,7)/11$ every coefficient is divisible by **40**. The largest forced divisors near $(1,2,4)/7$ and $(5,13,19)/37$ are **3** and **14**. In the first and third examples with divisors 40 and 14 the conclusion holds with every part proportion within $10^{-6}$ of its target. For divisor 3 at $(1,2,4)/7$, the tolerance is $1/100$. Any three rows of one weight in the same indicated neighborhood inherit the corresponding common divisor. The [equal-part theorem](equal_part_cones.md) gives a sharper classification on $(1,1,m)$, including neighborhoods forcing arbitrary prescribed powers of any prime.

## A finite graph read from the leading digits

Normalize an actual row and its weight by the same $p^n$, so $1/p\le N/p^n<1$. Pad with zeros indefinitely on the right. Read corresponding digits of the parts from most to least significant. Let $v_i$ be their prefix integers at a given depth, $S=\sum v_i$, and define discrepancies $z_i=Dv_i-s_iS$. The graph retains
$$
               s_i-D\le z_i\le2s_i,\qquad \sum z_i=0.             \tag{1}
$$
A state stores these discrepancies, the carry $c\in\{0,1,2\}$ that the next column must send left, and the remaining valuation budget. Start at zero discrepancies, carry zero, and budget $B$. For part digits $d_i$, put $t=\sum d_i$. Choose the incoming carry $c'\in\{0,1,2\}$ from the unread suffix. Require $0\le t+c'-pc<p$, strictly positive at the first column. Update by
$$
              z_i'=pz_i+Dd_i-s_it,\qquad c\longmapsto c',
$$
retain (1), and subtract $c$ from the budget without allowing it to become negative. Kummer's valuation counts the amounts of outgoing carry, including a carry of two twice. Equivalently the contribution of a layer is the factorial floor of the total minus the sum of the three part floors. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607)

For fixed $S\bmod D$, each of $z_1,z_2$ lies in one prescribed class modulo $D$ in an interval of length $D+s_i<2D$. There are at most four choices together, and $z_3$ is determined. This gives the claimed state bound. At budget zero the unread part tails sum to at most one, so the sharper upper bound $z_i\le s_i$ also holds in the exact limit. The same integer-rounding argument below gives it for the finite tolerance. The implementation chooses the digit of the total first. For fixed $t$, (1) leaves at most two choices for either of the first two part digits.

For an explicit tolerance, set $H=12D(B+1)+1$. Any row with
$$
          \max_i|x_i/N-s_i/D|\le 1/(D p^H)
$$
satisfies (1) throughout the first $H$ columns. Indeed, at depth $j\le H$ write the scaled parts $p^j x_i/p^n=v_i+f_i$ with $0\le f_i<1$. Then $z_i$ equals
$$
 D(p^jN/p^n)(x_i/N-s_i/D)+(s_i-D)f_i+s_i\sum_{l\ne i}f_l.
$$
The first term has absolute value less than one. The other terms are between the bounds (1), and $z_i$ is integral. If the valuation is at most $B$, the row therefore supplies $H$ states after the first edge, so some state repeats. The initial condition on the leading digit is retained when that cycle is repeated.

Conversely, if there is a reachable cycle after the leading digit, take a prefix to it and repeat it. Its cost is zero because the remaining budget returns to its old value. Bounded discrepancies ensure the infinite part words represent proportions exactly $s/D$. The digit equation, summed with place weights, shows their total equals the infinite total word, whose first digit is nonzero. For each truncation, replacing the prescribed carry from below by zero can only reduce the carries to its left. Its valuation is at most $B$. If the integer vector of prefix words is $a$ and that of one cycle is $b$, of length $h$, truncating after $n$ cycles gives
$$
        a+\big((p^h-1)a+b\big)(p^{hn}-1)/(p^h-1).
$$
The positive vector in parentheses is proportional to $s$, since the directions converge to $s/D$. Primitivity makes it an integer multiple $Ls$. This proves both alternatives and justifies deciding them by an ordinary directed-cycle search.

Only primes dividing $K(s)$ can enter $G(s)$, and their exponents are bounded by those in $K(s)$. The neighborhoods for the finitely many prime powers can be intersected. In fact $G(s)\mid F(s)/D$, where $F(s)=\gcd_{w\ge1}K(ws)$ is the [fixed divisor on the exact direction](fixed_direction_divisors.md). To see this, decrease one coordinate of $p^nws$ by one. Its lowest $n$ digits are all $p-1$ in that slot and zero in the others. The valuation is exactly that of $K(ws-e_i)$, with zero coordinates allowed in this last expression. The gcd over $i$ of $K(ws-e_i)=(s_i/D)K(ws)$ is $K(ws)/D$. Taking all $w$ proves the claim.

## What limiting proportions can decide at cofactor five

At the two [remaining odd prime powers below $2^{23}$](factor_five_exponents.md), a hypothetical counterexample requires a row approaching $s/(5Q)$ for $s=(3Q,z,2Q-z)$. For **each** $p=2,3,5,Q$ there are rows arbitrarily close to this direction whose coefficients avoid $p$. Thus no neighborhood alone can certify any of these four required primes here. The rows witnessing the different primes can have different weights, and the original arithmetic restriction on the scalar remains essential to the open problem.

| $Q$ | $z$ | Binary prefix, period lengths | Ternary | Quinary |
|---:|---:|---:|---:|---:|
| 1,467,779 | 2,918,400 | 23, 1 | 714, 19 | 321, 324 |
| 2,823,523 | 5,636,100 | 24, 1 | 971, 170 | 64, 160 |

For these positive witnesses the search tightens the first two upper bounds to $s_i-1$. Each resulting prefix and cycle is independently checked by proportionality of its integer slope vector and factorial floors on the first three repetitions. Both $Q$ are prime. The row $(3QT,zT-1,(2Q-z)T)$ at $T=Q^n$ avoids $Q$, since its upper parts $(3Q,z-1,2Q-z)$ add without a base-$Q$ carry.

## Reproduction

Run **python3 check_rational_cones.py** for 1,055 comparisons with the separate equal-part graph, checks of the three displayed exact divisors, and factorial-floor checks of the returned periodic rows. Running **python3 check_rational_cones.py 1 3 7** prints the exact divisor 40. To test one prime directly, pass the three slopes, the prime, and the maximum allowed valuation. For example, **python3 check_rational_cones.py 1 3 7 2 2** excludes valuations at most two. Changing the final argument to 3 returns a prefix and cycle attaining that bound. [General solver](check_rational_cones.py)

The larger practical tolerances are verified independently by retaining actual leading part integers. At a prefix with total digit word $n$ and part word $v_i$, a part share in $[l_i,u_i]$ requires $v_i<(n+1)u_i$ and $n l_i<v_i+1$. The total prefix equals the sum of the part prefixes plus the carry from below. Try all part digits and carries up to the specified budget and discard incompatible prefixes. The required trees are empty by depth 13 or earlier for the three examples above. [Explicit neighborhood checker](check_equal_part_cones.py)
