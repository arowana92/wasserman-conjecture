"""Verify all consecutive binomial cuts through a given quotient m=a(a+1)/N."""
import argparse
from math import isqrt
from pathlib import Path
from check_cross_cut_audit import common_primes, divisors


def scan(limit):
    count = largest = 0
    parity = set()
    for m in range(3,limit+1):
        for predecessor in divisors(m*(m-2)):
            n = predecessor+1
            if n <= m+1:
                continue
            root = isqrt(1+4*m*n)
            a = (root-1)//2
            if a+1 >= n or a*(a+1) != m*n:
                continue
            count += 1
            largest = max(largest,n)
            if (a & (n-a)) and ((a+1) & (n-a-1)):
                continue
            p = min(common_primes(n,a,a+1))
            parity.add((m,n,a,a+1))
            print(m,n,a,a+1,"least common endpoint prime",p)
    return count,largest,parity


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit",type=int,nargs="?",default=10**6)
    parser.add_argument("--audit-log",type=Path)
    args = parser.parse_args()
    count,largest,parity = scan(args.limit)
    print("Systems",count,"parity",len(parity),"largest",largest)
    if args.limit == 10**6:
        assert (count,len(parity),largest) == (1813,6,9896020401)
    if args.audit_log:
        *records, summary = args.audit_log.read_text().splitlines()
        assert summary == f"Systems {count} parity {len(parity)} largest {largest}"
        assert len(records) == len(parity)
        assert {tuple(map(int,row.split())) for row in records} == parity
        print("All parity survivors agree with the independent q-first enumeration.")
