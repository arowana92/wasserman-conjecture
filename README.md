# Wasserman’s common-divisor conjecture for trinomial coefficients

**Work in progress.** This repository studies the open general case and records partial proofs and reproducible computations.

For an integer $N\ge3$, take any three triples of positive integers summing to $N$. Wasserman conjectured that a prime divides all three numbers $N!/(a!b!c!)$ obtained from the triples.

The conjecture has been checked here for every $3\le N\le5000$. At each weight the table gives a smallest possible set of candidate common primes. The only weights in the checked range requiring five primes are 699 and 717. A separate proof gives an infinite family of composite weights with certificates of size at most seven.

- [Proof for all sufficiently large odd weights with three binary ones](three_binary_digits.md)
- [The conjecture throughout growing intervals above every power of two](power_intervals.md)
- [Proof for every weight $5\cdot2^k+1$ and eventual results for $m2^k+1$](binary_weights.md)
- [Research notes, finite theorem, and reproduction commands](RESEARCH.md)
- [Reduction and verification through one million for selections containing $(1,1,N-2)$](endpoint_reduction.md)
- [Endpoint part-size bound and an exclusion at 2,099,636,978 weights](endpoint_factors.md)
- [Restrictions on counterexamples from weighted factorial products](weighted_obstruction.md)
- [Proof for an infinite family](bounded_certificates.md)
- [Prime certificates through 5,000](certificates_5000.txt)
- [George Bergman’s paper and Proposition 15](https://arxiv.org/abs/0806.0607)
