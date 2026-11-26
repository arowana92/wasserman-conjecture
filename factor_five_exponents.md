# All binary exponents at bounded odd prime powers

Let $Q<2^{23}$ be an odd prime power and $k\ge1$. The trinomial common-divisor conjecture holds at $N=5Q2^k+1$ except possibly on the following two progressions, where $k=k_*+st$ and $t\ge0$. The parameter $z$ determines the only possible row below.

| $Q$ | $z$ | $k_*$ | $s$ |
|---:|---:|---:|---:|
| 1,467,779 | 2,918,400 | 2,337,019 | 2,935,556 |
| 2,823,523 | 5,636,100 | 3,948,614 | 5,647,044 |

On both progressions the conjecture holds for all $t<3^{256}$, covering every $k\le10^{128}$. In any interval of $H\ge1$ consecutive exponents, at most $O(H^{26/625})$ values of $k$ can give exceptions on the two progressions combined. On the second progression alone the bound is $O(H^{1/25})$. In particular the conjecture holds for every positive exponent at every odd prime power $Q<1{,}467{,}779$.

## Determining the row from a binary subset

Put $P=2^k$. The [factor-five reduction](predecessor_five_seven.md#full-factor-five) forces $Q=\ell^e$ with $\ell>5$ prime, $Q\equiv3\bmod8$, at least seven binary ones in $Q$, and no common bit of $3Q$ and $2Q$. The endpoint row occurs. One more row has three equal parts $N/3$, with $3\parallel N$. The last row
$$A=(3QP,zP,(2Q-z)P+1)$$
must avoid divisibility of its coefficient by $2,3,5,\ell$. The bits of $z$ are a proper nonempty subset of the bits of $2Q$, $15\mid z$, and $zP\equiv1\bmod Q$. Also $PQ\equiv7\bmod15$, so $k$ has a prescribed residue $k_0\bmod4$.

For any candidate $z$ put $w=z^{-1}\bmod Q$, taking the least nonnegative residue. For each $g=\ell^j$, $1\le j\le e$, the two groups $3QP,2QP+1$ must add without a carry at modulus $Qg$. Their least residues there are $Q[3w]_g$ and $Q[2w]_g+1$, where $[x]_g$ is the least residue. Thus $[3w]_g+[2w]_g<g$, equivalently $[3w]_g\le[5w]_g$. Next solve $16^u=w2^{-k_0}\bmod Q$. The unit group is cyclic, so exponentiation by $\operatorname{ord}_Q(16)$ tests membership in $\langle16\rangle$. A baby-step giant-step search gives $u$ when it exists. The resulting $k$ form one progression with stride $s=4\operatorname{ord}_Q(16)$.

## Testing every higher digit

For each progression and each $p=3,5$, start with modulus $m=p^{v_p(2^s-1)}$, increasing it to 9 at $p=3$ if needed. Enumerate all possible residues $W=P\bmod m$ by iterating multiplication by $2^s\bmod m$. Discard $9\mid N$ at $p=3$, and test the low digits of the three parts by Kummer. Write $(a_1,a_2,a_3)=(3Q,z,2Q-z)$ and $(b_1,b_2,b_3)=(0,0,1)$. The starting state has entries $S_i=\lfloor(a_iW+b_i)/m\rfloor$.

For each possible next input digit $v$ of the scalar, the three output digits are $[S_i+a_iv]_p$. If their sum is below $p$, take the transition to the state with entries $\lfloor(S_i+a_iv)/p\rfloor$. The states obey $0\le S_i\le a_i$ by induction. A finite scalar can end at any state whose three entries add without further carries, by appending zeros. Therefore an exhaustive search from the starting states that fails to reach any such state rules out the entire progression. Allowing every possible higher scalar digit preserves every actual power of two.

Two programs enumerate the progressions independently. The Python program generates binary strings and subsets, solves the modular logarithms and tracks all three carries. The C++ program sieves the integers, iterates inverse powers until they repeat, and tracks the first two parts and the total. In each column those two digits must sum to at most the digit of the total. Their counts agree at both cutoffs.

| Strict upper bound for $Q$ | Eligible prime powers | Progressions | Excluded at 3 | Further exclusions at 5 | Remaining |
|---|---:|---:|---:|---:|---:|
| $2^{21}$ | 122 | 1,194 | 1,173 | 20 | 1 |
| $2^{23}$ | 269 | 3,912 | 3,860 | 50 | 2 |

## The remaining progressions

At the first exponent $k_*$ on the first progression, the three residues at modulus 27 are $(3,21,18)$ and the total has residue 15. At modulus 125 the residues are $(56,75,5)$ and the total has residue 11. On the second progression the starting exponent gives $(3,21,9)$ against 6 at modulus 27, and $(21,0,15)$ against 11 at modulus 25. Each overflow excludes the corresponding exponent.

A productive state is one from which a terminating state can be reached. The ternary graphs have the following sizes. The outdegrees count transitions to productive states, including loops.

| $Q$ | Reachable | Terminating | Productive | Outdegree one | Outdegree two |
|---:|---:|---:|---:|---:|---:|
| 1,467,779 | 151,397 | 21 | 29,656 | 28,350 | 1,306 |
| 2,823,523 | 199,807 | 46 | 36,778 | 35,089 | 1,689 |

Starting at each productive state with one empty path, the recurrence summing counts at its productive successors gives the respective maximum counts of length-2,048 paths
$$M_1=40879096506588498826141913655374463598481,\qquad
  M_2=878909565367958979420167123851270494639.$$
Integer comparison gives $M_1^{625}<3^{26\cdot2048}$ and $M_2^{25}<3^{2048}$. Every admissible prefix comes from a path through productive states. Because each productive state has an outgoing edge in this set, every shorter path extends to length 2,048. For a fixed starting state and a fixed shorter length, the number of paths is therefore bounded by the corresponding $M$. Splitting the input into blocks, the number of admissible length-$j$ prefixes is $O(M^{\lceil j/2048\rceil})$, including the finitely many starting residues.

Write $v=v_3(2^s-1)$, which is 1 on the first progression and 2 on the second. Lifting gives $v_3(2^{su}-1)=v+v_3(u)$ for $u>0$. Thus the residue of $P\bmod 3^{j+v}$ determines $t\bmod 3^j$. Every interval of at most $3^j$ successive $t$ has pairwise distinct residues. The first two ternary digits are handled by the starting states. Take $j=\lceil\log_3 H\rceil$ and apply the path bound to obtain the stated estimates, uniformly in the location of the interval. An interval of $H$ consecutive $k$ meets at most $H$ consecutive values of $t$ on each progression.

For the explicit initial range, lift one ternary digit of $t$ at a time. Extending $t$ by $u3^j$ multiplies the scalar residue by $(2^{s3^j})^u$. Use arithmetic modulo $3^{256+v}$ and test the three least residues at layer $3^{j+1+v}$ after each extension. Respectively 83,502 and 205,161 values $0\le t<3^{256}$ satisfy $9\nmid N$ and pass every ternary layer through $3^{256+v}$. Independently, counting paths of $256+v-2$ digits from the starting states, including all reachable states, gives the same counts. The checker recomputes $2^k\bmod 3^{256+v}$ by direct modular exponentiation for every reported survivor and checks all its ternary digits.

For each survivor compute $2^k\bmod 5^{10}$ directly. Compare the floor of the total divided by each $5^i$ with the sum of the three floors. The first positive difference occurs at the following layer in every case.

| $Q$ | $i=2$ | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1,467,779 | 33,467 | 46,664 | 2,073 | 1,035 | 213 | 13 | 35 | 0 | 2 |
| 2,823,523 | 81,956 | 82,110 | 26,400 | 9,162 | 4,445 | 881 | 160 | 34 | 13 |

## Reproduction

Run the following with C++17, Python 3 and python-flint. The Python programs use assertions to compare exact counts.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_factor_five_rays.cpp -o .work/check_factor_five_rays
.work/check_factor_five_rays 23
python3 check_factor_five_rays.py 23
python3 check_factor_five_tail.py
~~~
