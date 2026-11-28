# The 4096th power with an independent binary exponent

Write $P=65537$, $R=2^b$, $T=q^{4096j}$, and $N=1+PRT$, with $q$ prime and $b,j\ge1$. The three-row conjecture holds for every such input with **$b<72$**. In the unrestricted range it holds unless $4\mid b$ and the following two congruences both hold:
$$
 2^{b+14}T\equiv-1\bmod P,\qquad
 b\equiv j\,Q_P(q)+c\bmod P,\quad c\in\{-14,2,18,34,50\}, \tag{1}
$$
where $Q_P(a)=(a^{P-1}-1)/P\bmod P$ for every $a$ coprime to $P$. For fixed $q,j$ at most five classes of $b$ modulo $32P=2,097,184$ survive (1). The first congruence either has no solutions or specifies a single class modulo 32. Coincident prime bases are already covered by the [two-prime predecessor theorem](all_weights_endpoint_reduction.md#consequences).

## Two consequences of the digit at position one in base $P$

At a failure the [half-exponent reduction](independent_fermat_exponents.md#the-remaining-conditions-at-half-the-exponent) supplies $4\mid b$, an odd $j$, a unitary scale $U=13$ and numerator $x=49153=P-2^{14}$. In the row which avoids $2,q,P,13$, the part with residue one modulo $P$ is $xRT$. Its part of residue one modulo $R$ is $PzT$ for some positive $z$, divisible by 13. Put $S=2^{b+14}$. The three actual parts are
$$
               xRT,\quad PzT,\quad (S-Pz)T+1. \tag{2}
$$
The third part is a positive multiple of $PR$. Indeed the distinct full powers $P,R,T$ have distinct exceptional positions, each of residue one. In particular $0<Pz<S$.

Since $xRT\equiv1\bmod P$ and $4x=3P+1$, we have $RT\equiv4\bmod P$. The units of the total in base $P$ are one, and its next digit is four. The digit in this second position of $PzT$ is $zT\bmod P$. Hence
$$
            zT\bmod P\in\{0,1,2,3,4\},\qquad
            4z\equiv rR\bmod P\quad(0\le r\le4). \tag{3}
$$

The digit in the same position of $xRT$ is also at most four. Write $RT=4+Pu$, so $(xRT-1)/P=3+xu$. The five possibilities for this digit give $u\equiv4s\bmod P$, $-3\le s\le1$. It follows that
$$
                        RT/4\equiv1+Ps\bmod P^2
$$
using the inverse of four modulo $P^2$. The Fermat quotient adds under multiplication, as follows by multiplying $a^{P-1}\equiv1+P Q_P(a)\bmod P^2$. Also $Q_P(2)=-4096\bmod P$, by expanding $2^{65536}=(1-P)^{4096}$. Taking the Fermat quotient on both sides yields
$$
     -4096(b-2)+4096j Q_P(q)\equiv -s\bmod P.
$$
Since $4096^{-1}\equiv-16\bmod P$, this proves (1).

## A finite contradiction independent of the prime $q$

The binomial cut with first part $xRT$ is odd exactly when $xT$ and $2^{14}T$ have disjoint bits. The other two parts must also add freely to $1+ST$. The third is even, so below the power $S$ the entire unit residue belongs to $PzT$. In particular
$$
       PzT\equiv1\bmod S,\qquad B=(PzT-1)/S
       \text{ has all its occupied bits among those of } T . \tag{4}
$$
This follows by dividing the two parts by $S$ and discarding their units residue, leaving a carry-free sum $B+(T-B)=T$. The preceding reduction gives $T\equiv1+2^{14}\bmod2^{15}$, so $z\equiv1+2^{14}\bmod2^{15}$ by (4).

For each multiple of four through 68, enumerate all $0<z<S/P$ with these low digits, $13\mid z$, and (3). Its forced input prefix is $t=(Pz)^{-1}\bmod S$. Retain it if the bottom $b+14$ bits of $xt$ and $2^{14}t$ are disjoint. This gives the second column below.

| $b$ | Candidate parts $z$ | Maximum additional digits until all extensions fail |
|---:|---:|---:|
| 4 through 52, in steps of four | 0 | 0 |
| 56 | 2 | 5 |
| 60 | 16 | 8 |
| 64 | 74 | 40 |
| 68 | 392 | 26 |

To extend a prefix by a bit at $S2^i$, try both choices. Its new digit must preserve disjointness of $xT$ and $2^{14}T$ at that position. It also fixes bit $i$ of $(PzT-1)/S$, which may be occupied only if bit $i$ of $T$ is occupied. All branches die at the depths in the table. The search includes zero padding of finite positive inputs, so the obstruction holds for all $q,j$.

Two implementations enumerate the starting parts. One generates the allowed low prefixes of $T$ and inverts them. The other branches directly on the low bits of $z$ and lifts its inverse using Newton's identity. To check the extensions independently, enumerate digits of $B$ instead. Its low $i$ digits specify $T=(1+SB)(Pz)^{-1}\bmod S2^i$. At each new modulus, compare the sum of integer quotients of the three actual parts in (2) with the quotient of their total. Do the same for the independent binary cut $(xT,2^{14}T)$. These are Legendre's zero-defect tests at that layer. Both checks give the same depth at every individual starting part.

From the repository root, run the commands below to compare the two complete sets of parts and repeat both extension proofs.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 verify_fermat_half_seeds.cpp -o .work/verify_fermat_half_seeds
.work/verify_fermat_half_seeds | python3 check_fermat_half_prefixes.py --compare-seeds
~~~

The finite range establishes a theorem at unbounded weights. The five-class condition still permits larger binary exponents.
