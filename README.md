# Wasserman’s common-divisor conjecture for trinomial coefficients

**Work in progress.** This repository studies the open general case and records partial proofs and reproducible computations.

For an integer $N\ge3$, take any three triples of positive integers summing to $N$. Wasserman conjectured that a prime divides all three numbers $N!/(a!b!c!)$ obtained from the triples.

Every counterexample must contain $(1,1,N-2)$ up to permutation. This reduces the full conjecture at every weight to proving that any two coefficients share a prime factor of $N(N-1)$. In particular the conjecture holds whenever $N-1$ has at most two distinct prime divisors.

The conjecture holds at a set of weights of natural density one. The number of possible exceptions up to $X$ is $O(X^{15/16})$.

The full assertion is proved here for every $N=5\cdot2^k+1$, for every weight with base-$p$ digit sum at most five for some prime $p\mid N-1$, and throughout growing intervals above powers of two.

The conjecture has been checked here for every $3\le N\le1{,}000{,}000{,}000$. Through 5,000 the table gives a smallest possible set of candidate common primes. The only weights through 5,000 requiring five primes are 699 and 717. For all sufficiently large $X$, at least $cX/(\log X)^2$ pairs of consecutive weights in $[X,2X]$ have minimum certificate size exactly four at both weights.

- [Exponent spacing for two digit blocks and an exceptional-set bound for three free powers](two_block_spacing.md)
- [Three independently placed prime powers, with complete results for five small offsets](two_high_powers.md)
- [Density one of two-prime-power cofactors at every exponent, and an exhaustive binary check through factor 1,000](binary_two_factors.md)
- [A sharp count for binary carry-free affine parts and bounds at every prime-power cofactor](binary_affine_counts.md)
- [Two consecutive prime powers and an affine third factor, with 8,192 families checked at every input](consecutive_power_partners.md)
- [Reciprocity between predecessor primes and a search of prescribed CRT rows](crt_reciprocity.md)
- [A digit test for three predecessor primes and two polynomial prime-power factors](cofactor_digit_gap.md)
- [Polynomial weights at almost all exponents of every prime, with power-saving bounds for prime bases](polynomial_row_divisors.md)
- [Products from growing blocks of consecutive primes, their successors and an explicit gap bound](translated_products.md)
- [Prime bases that work at every exponent and cofactor through 1,000](prime_base_density.md)
- [Shifted row denominators and all weights $d p^{2e}+1$ through cofactor 1,000](shifted_directions.md)
- [Polynomial weights with a common coefficient factor](polynomial_content.md)
- [Polynomial weights at arbitrary prime powers, with effective resultant bounds](polynomial_prime_powers.md)
- [Whole intervals on generic polynomial families and 1,021 consecutive weights at every binary exponent](polynomial_intervals.md)
- [Quadratic prime-power families with a discriminant criterion and complete checks for coefficients through 15](quadratic_powers.md)
- [Unrestricted prime-power exponents through offset 191, with four possible rays through 256](successor_exponents.md)
- [Two consecutive weights near a product of prime powers, for all offsets through 1,000](successor_powers.md)
- [Common-divisor bounds from determinants and density in short intervals after every fixed binary prefix](projected_rows.md)
- [Classification of shifted proportional rows through binary prefix 4,095](shifted_binary_prefixes.md)
- [All weights $p^a(p^j+d)$ for $j\ge a\ge1$, $1\le d\le1000$, and $p\nmid d$](close_prime_powers.md)
- [Optimal four-prime certificates at consecutive weights](four_prime_pairs.md)
- [A spacing theorem for six base-$p$ digit units, with finite certificates](six_digit_gaps.md)
- [Proof of the reduction to the endpoint at every weight](all_weights_endpoint_reduction.md)
- [Simultaneous carry tests through a thousand-digit bound and counts at unbounded weights](carry_intersections.md)
- [Two predecessor prime powers and the family $6p^a q^b+1$](two_prime_powers.md)
- [The full factors five and seven, with bounds for two varying parameters](predecessor_five_seven.md)
- [Three prime factors of the predecessor with a full power at most four](three_factors_small_power.md)
- [Four nearby prime multiples and restrictions on parts 1 and 3](near_row.md)
- [Certified intervals at offsets from powers of a prime](prime_power_offset_intervals.md)
- [Quantitative endpoint gcd bound from two parts of a row](endpoint_pair_bound.md)
- [Endpoint reduction through $10^{20}$ and a critical-exponent bound for small parts](critical_small_parts.md)
- [Complete finite proof through one billion and the prime-power layer bound](full_range_1b.md)
- [Uniform theorem for weights with two separated binary blocks](separated_binary_blocks.md)
- [Polynomial families evaluated at powers of two](polynomial_binary_families.md)
- [Density one along every fixed odd-cofactor prime-power ray](prime_power_ray_density.md)
- [Reduction to the endpoint for every bounded binary digit sum](bounded_binary_weight.md)
- [Exceptional exponent patterns among odd weights with a fixed number of binary ones](odd_sparse_endpoint.md)
- [Infinite families from a large prime power dividing the predecessor](prime_power_cofactors.md)
- [Density-one prime lists with mean size below 35.32](typical_certificates.md)
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
