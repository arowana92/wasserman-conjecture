# Common exponents for growing sets of primes and cofactors

Let $H$ tend to infinity, and take **any** increasing arithmetic progression of $H$ positive exponents. Outside at most
$$
 O\!\left(H\exp\big(-\sqrt{\log H}/6\big)\right)
$$
terms $e$, Wasserman's full conjecture holds simultaneously at every weight
$$
 1+dq^e,\qquad q\le H^{0.44}\ \text{prime},\quad
 1\le d\le \exp\big(\tfrac12\sqrt{\log H}\big).
$$
More generally the bases can be 2 and **any prescribed collection of at most $H^{0.44}$ odd primes**, of unrestricted size. The constant is absolute. The progression therefore contains a run of $\gg\exp(\sqrt{\log H}/6)$ consecutive terms all satisfying the simultaneous assertion, by dividing the good terms into gaps between bad ones.

## The base two

Write $D=\lfloor\exp(\tfrac12\sqrt{\log H})\rfloor$. At an exceptional weight $1+d2^e$ with $d\le D$, the [endpoint assignment](small_weight_transfer.md#specializing-the-two-rows) supplies a row
$$
 (a2^e,b2^e,c2^e+1),\qquad a+b+c=d,\quad a,b,c\ge1 .
$$
Its coefficient is odd, and it avoids every prime factor of a divisor $U\ge3$ of the weight, with $U\mid a,b$. In particular it avoids some odd prime $r\mid\gcd(a,b)$, with $r\nmid d,c$. The binary digits of $a,b,c$ are disjoint since the row avoids 2. Over *all* $d\le D$ there are at most $4^{\lceil\log_2(D+1)\rceil}=O(D^2)$ such ordered triples. Each bit belongs to one of the three slopes or to none. There are $O(\log(2D))$ candidate primes per triple.

Fix a triple and prime $r$. As functions of the input $Y=2^e$, the three parts and their total have slopes $a,b,c,d$. Exactly two are units modulo $r$. Let $J=\max v_r(a_i)$ over these four slopes, so $J\le\log D/\log r$. From any fixed compatible low input prefix the number of admissible strings of $l$ subsequent digits is at most
$$
 e^J B_r^l,\qquad
 B_r=\begin{cases}
 (2r+1)/3 & r<101,\\
 (27/50)r & r\ge101 .
 \end{cases} \tag{1}
$$
Here $e^J$ denotes the exponential function of $J$. To check the constant, use the [odd-base affine digit proof](binary_affine_counts.md#an-odd-prime-base). With exactly two unit slopes and the other two delayed, the sum of relative prefix weights at each digit is at most
$$ \frac{2\exp(1)}{3\sqrt3(1-\exp(-2/r))}\le \frac{\exp(1)(r+2)}{3\sqrt3} $$
for $r\ge5$, using the pending digit weight $\exp(-b/r)$ for pending digit sum $b\le (r-1)J$. This gives the first rate at $5\le r<101$, as in the linked proof, and the second at $r\ge101$. For $r=3$ use weight $(5/4)^{-b}$, whose terminal loss is at most $(5/4)^{2J}\le\exp(J)$. After stripping a fixed low input prefix, divide the four affine expressions by its place value, rounding down their offsets. At a compatible prefix no carry leaves the omitted columns, so the new offsets for the parts still sum to that for the total. The slopes, delays and bounds are unchanged. This proves (1) for arbitrary starting prefixes.

Split the exponent indices into $r-1$ classes. In each, write $Y_i=2^{A+Si}$ with $r-1\mid S$. For $i\ne j$,
$$ v_r(Y_i-Y_j)=h+v_r(i-j),\qquad h=v_r(2^S-1), $$
by [lifting the exponent](moving_power_bases.md#reading-the-digits-that-actually-vary). All inputs in the class have a common prefix of length $h$. If the class has $M$ terms, take $l=\lfloor\log_r M\rfloor$. Any particular next $l$ digits occur at most $\lceil M/r^l\rceil$ times, since they fix the index modulo $r^l$. With $s_r=\log_r B_r$, (1) bounds the number of surviving indices by
$$ 2\exp(J) M^{s_r}r^{1-s_r}. $$
Summing the classes costs at most another $(r-1)^{1-s_r}$ by concavity. Both losses are bounded absolutely, since $r^{1-s_r}=r/B_r\le 50/27$. Thus the count for this triple and prime is $O(\exp(J)H^{s_r})$, independent of the initial exponent and progression step.

For the finitely many odd primes $r<101$, $s_r$ is bounded strictly below one. Their total, even after counting all the triples and delay factors, is $O(H^{1-\eta})$ for some fixed $\eta>0$. For $r\ge101$ we have $27/50<\exp(-3/5)$, hence
$$
 \exp(J)H^{s_r}
 \le H\exp\!\left(-\frac{(3/5)\log H-\log D}{\log r}\right)
 \ll H\exp\big(-(6/5)\sqrt{\log H}\big).
$$
Here $r\le D$ and the numerator is positive for large $H$. The $O(D^2\log(2D))$ total candidates cost at most $\exp(\sqrt{\log H}+O(\log\log H))$. This proves the exceptional bound at base 2.

## The odd bases

For each odd prime $q$, the [uniform theorem for moving powers](moving_power_bases.md) bounds the number of terms that fail for *any* $d\le D$ by $O(D^{4-\alpha}H^\alpha)$ with an absolute constant and $\alpha<0.552$. Summing over the prescribed bases gives $O(H^{0.992+o(1)})$, which is absorbed by the stated error. The proof handles progression steps divisible by arbitrarily large prime powers in both cases.

Run `python3 check_simultaneous_exponents.py` for exact slope counts, branching constants, and residue tests across shifted exponent progressions.
