"""Short integer certificates for complete carry-free extensions in every modulus."""

from math import gcd


# Each entry gives the radix, its initial modulus and residue, an accepting
# suffix, and pairs of return words (number of digits, first word, second word).
# A word is written as its integer value, read from the units place upwards.
CERTIFICATES = {
    (12291, 384, 2749): (
        (2, 64, 43, 50, (
            (23, 5636146, 5636298),
            (43, 5909981265970, 5910087532594))),
        (3, 729, 200, 150599, (
            (65, 2826377662422966618940684896161, 2826377662422966625791807405131),
            (66, 8479132987268899877158618788962, 8479132987268899877494706013962),
            (68, 76312196885419734646722690030440, 76312196885420098711392830106206),
            (68, 76312196885419734646722690030440, 76312196885420098899443848301915))),
        (241, 241, 150, 172, ((3, 8712322, 8712386), (3, 8712322, 8712467)))),
    (11415, 320, 12617): (
        (2, 32, 7, 21312, (
            (61, 504426042086150976, 504448925906785088),
            (77, 33059564800323165901632, 33059564808228464251712))),
        (5, 25, 22, 668647716, (
            (21, 419619800586231, 419620044961231),
            (22, 2098100212286671, 2098100212317646))),
        (761, 761, 283, 0, ((2, 215363, 215368), (2, 215363, 215370)))),
    (15587, 8448, 3869): (
        (2, 256, 203, 9330507358952, (
            (136, 69077408609100948256931573823019494384360,
             69077408609100948257551458820839797662440),
            (143, 8841908301964921376966625824549043717776104,
             8842268894136213485992887323759203428049640))),
        (11, 121, 114, 3658516, (
            (11, 268803071387, 268803322217),
            (12, 2956813398400, 2956836551802),
            (12, 2956813398400, 2956836604437))),
        (109, 109, 2, 17, ((3, 23779, 23793), (4, 2590402, 2590729)))),
}


def layer_check(entries, base, through):
    layer = base
    while layer <= through:
        assert sum(entries) // layer == sum(x // layer for x in entries)
        layer *= base


def verify():
    pairs = 0
    for slopes, rows in CERTIFICATES.items():
        for p, power, low, suffix, words in rows:
            def entries(t):
                return (slopes[0]*t, slopes[1]*t, slopes[2]*t+1)

            assert 0 < low < power and gcd(low, p) == 1
            part = power
            while part % p == 0:
                part //= p
            assert part == 1
            finished = entries(low + power*suffix)
            layer_check(finished, p, sum(finished))
            state = [x // power for x in entries(low)]
            divisor = 0
            for length, first, second in words:
                block = p**length
                assert 0 < first < second < block
                for v in (first, second):
                    row = entries(low + power*v)
                    layer_check(row, p, power*block)
                    assert [x // (power*block) for x in row] == state
                    terminal = entries(low + power*(v + block*suffix))
                    layer_check(terminal, p, sum(terminal))
                divisor = gcd(divisor, second-first)
                pairs += 1
            while divisor % p == 0:
                divisor //= p
            assert divisor == 1
        w, P, u = rows[0][1], rows[2][0], rows[1][0]
        assert sum(slopes) == w*P
        assert gcd(slopes[0], slopes[1]) == u
        assert slopes[0]*rows[0][2] % w == 1 and slopes[1]*rows[2][2] % P == 1
        assert (1 + w*P*rows[1][2]) % u == 0
        assert (1 + w*P*rows[1][2]) % (u*u) != 0
    print(pairs, "pairs of return words verified across nine radix and cofactor choices")


if __name__ == "__main__":
    verify()
