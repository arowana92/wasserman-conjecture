# Wasserman’s common-divisor conjecture for trinomial coefficients

**Work in progress.** This repository studies the open general case and records partial proofs and reproducible computations.

For an integer $N\ge3$, take any three triples of positive integers summing to $N$. Wasserman conjectured that a prime divides all three numbers $N!/(a!b!c!)$ obtained from the triples.

The conjecture holds at a set of weights of natural density one. The number of possible exceptions up to $X$ is $O(X^{15/16})$.

The full assertion is proved here for every $N=5\cdot2^k+1$, for all sufficiently large odd weights with three binary ones, and throughout growing intervals above powers of two.

The conjecture has been checked here for every $3\le N\le1{,}000{,}000$. Through 5,000 the table gives a smallest possible set of candidate common primes. The only weights through 5,000 requiring five primes are 699 and 717. A separate proof gives an infinite family of composite weights with certificates of size at most seven.

- [Stronger exceptional-set bound and its factorization argument](endpoint_count.md)
- [Proof of density one with a power-saving exceptional set](almost_all_weights.md)
- [Complete finite proof for every weight through one million](full_range_1m.md)
- [All three-bit weights below $2^{128}$ and an endpoint reduction for larger weights](even_three_binary_digits.md)
- [Proof for all sufficiently large odd weights with three binary ones](three_binary_digits.md)
- [The conjecture throughout growing intervals above every power of two](power_intervals.md)
- [Proof for every weight $5\cdot2^k+1$ and eventual results for $m2^k+1$](binary_weights.md)
- [Research notes, finite theorem, and reproduction commands](RESEARCH.md)
- [Reduction and verification through one million for selections containing $(1,1,N-2)$](endpoint_reduction.md)
- [An odd-divisor theorem at all even weights below $2^{64}$](even_endpoint.md)
- [Endpoint part-size bound and an exclusion at 2,099,636,978 weights](endpoint_factors.md)
- [A lower bound on products of parts in different rows](parts_in_different_rows.md)
- [Restrictions on counterexamples from weighted factorial products](weighted_obstruction.md)
- [Proof for an infinite family](bounded_certificates.md)
- [Prime certificates through 5,000](certificates_5000.txt)
- [George Bergman’s paper and Proposition 15](https://arxiv.org/abs/0806.0607)
