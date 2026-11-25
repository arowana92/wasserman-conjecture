# Three prime factors of the predecessor, including a power at most four

**Theorem.** If $N-1$ has exactly three distinct prime factors and the smallest of its three full prime powers is at most 4, the full trinomial conjecture holds at $N$. In particular this includes every $N=4p^a q^b+1$ with distinct odd primes $p,q$ and $a,b\ge1$.

By the [all-weight reduction](all_weights_endpoint_reduction.md), it suffices to exclude the endpoint $E$. Write $N-1=w_1w_2w_3$ with the $w_i$ full powers of distinct primes. The [three-factor CRT argument](endpoint_factors.md) supplies a row $A$ avoiding all three bases. Each part has the form $e_i=((N-1)/w_i)k_i$ with $1\le k_i<w_i$, $e_i\equiv1\bmod w_i$. The row is assigned a unitary divisor $U$ of $N$ (that is, $\gcd(U,N/U)=1$), with $U\ge3$, whose prime bases it must also avoid. Since $U$ is coprime to $N-1$ and divides all parts, $U\mid k_i$ for each $i$. This excludes $\min w_i\le3$.

Suppose $w_i=4$ and write $t=(N-1)/4$. Then necessarily $U=k_i=3$. The row $A$ has a part $3t$, and its other two parts have sum $t+1$. All its parts are divisible by 3, so $3\mid t+1$. Put $k=v_3(t+1)\ge1$. Modulo $3^{k+1}$, the least nonnegative residue of $3t=3(t+1)-3$ is $3^{k+1}-3$, and that of $t+1$ is at least $3^k\ge3$. Their sum reaches the modulus. Adding the two summands $3t$ and $t+1$ therefore produces a carry in base 3. This forces 3 to divide $K_A=\binom{N}{3t}\binom{t+1}{e_j}$, contrary to the assignment of 3 to $A$.

More generally, whenever some $w_i=u+1$ with $u=p^e$ a prime power, the simultaneous assignment $k_i=U=u$ is impossible. Indeed $u\mid t+1$ for $t=(N-1)/w_i$, so $k=v_p(t+1)\ge e$. Modulo $p^{k+e}$, the residue of $ut$ is $p^{k+e}-u$, while the residue of its complement $t+1$ is at least $p^k\ge u$. This forces a carry in base $p$. The conclusion also holds when $k_i=u$ and the assigned $U$ is a positive power of the same $p$ provided $p^{v_p(t+1)}\ge u$.

Run `python3 check_three_factors.py` to test the carry obstruction on 808 exact binomial coefficients and reproduce the small CRT enumeration.
