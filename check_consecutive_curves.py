"""Match integral models against Cremona's published integral-point tables."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from math import isqrt
from pathlib import Path
from urllib.request import urlopen
from check_cross_cut_audit import common_primes

REVISION = "6fedd3ab463a70e52061ebc0b04832ad19fd2e01"
SOURCE = f"https://raw.githubusercontent.com/JohnCremona/ecdata/{REVISION}/intpts/"


def model(q):
    d, t = q, 1
    p = 2
    while p*p <= d:
        while d%(p*p) == 0:
            d //= p*p
            t *= p
        p += 1
    # When d is even and 8 divides q, the original x=d*m and y=d^2*t*a
    # are divisible by 4 and 8 since q divides m*(m-2), forcing m even.
    scale = 2 if d%2 == 0 and q%8 == 0 else 1
    a2, a3, a4 = -2*d//scale**2, d*d*t//scale**3, q*d*d//scale**4
    r = (1-a2)//3
    shift = -(a3//2)
    coefficients = (0, 3*r+a2, a3+2*shift, 3*r*r+2*a2*r+a4,
                    r**3+a2*r*r+a4*r-shift*shift-a3*shift)
    return d, r, scale, "["+",".join(map(str, coefficients))+"]"


def load_files(directory, download):
    paths = [directory/f"intpts.{start:05d}-{start+9999:05d}"
             for start in range(0,500000,10000)]
    if download:
        directory.mkdir(parents=True, exist_ok=True)
        def fetch(path):
            if not path.exists():
                with urlopen(SOURCE+path.name, timeout=60) as response:
                    data = response.read()
                tmp = path.with_suffix(".partial-"+path.suffix[1:])
                tmp.write_bytes(data)
                tmp.replace(path)
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(fetch,paths))
    return paths


def check(limit, paths):
    wanted = {}
    for q in range(1,limit+1):
        d, r, scale, equation = model(q)
        assert equation not in wanted
        wanted[equation] = (q,d,r,scale)
    matched, cases = {}, []
    for path in paths:
        with path.open() as rows:
            for row in rows:
                label, equation, coordinates = row.split()
                if equation not in wanted:
                    continue
                q, d, r, scale = wanted[equation]
                assert q not in matched
                matched[q] = label
                for coord in coordinates.strip("[]").split(","):
                    if not coord:
                        continue
                    m, rem = divmod(scale**2*(int(coord)+r), d)
                    if rem or m < 3 or m*(m-2)%q:
                        continue
                    n = 1+m*(m-2)//q
                    root = isqrt(1+4*m*n)
                    a = (root-1)//2
                    if a < 2 or a+1 >= n or a*(a+1) != m*n:
                        continue
                    witness = min(common_primes(n,a,a+1))
                    cases.append((q,m,n,a,witness))
    first_missing = next((q for q in range(1,limit+1) if q not in matched), limit+1)
    print("Matched curves:",len(matched),"first missing q:",first_missing)
    print("Matched q:",*sorted(matched))
    print("Admissible arithmetic cases:",len(cases),"largest m:",max(r[1] for r in cases))
    print("q=1:",matched[1],[r for r in cases if r[0] == 1])
    return matched,cases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory",type=Path,default=Path(".work/cremona_intpts"))
    parser.add_argument("--limit",type=int,default=1000)
    parser.add_argument("--download",action="store_true")
    args = parser.parse_args()
    found, cases = check(args.limit,load_files(args.directory,args.download))
    if args.limit == 1000:
        assert sorted(found) == [1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,25,27,
                                 32,36,40,56,64,72,128,384]
        assert len(cases) == 51 and max(row[1] for row in cases) == 56917632
