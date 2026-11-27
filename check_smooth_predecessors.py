"""Generate the complete divisor family and run the C++ predecessor sieve."""

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from math import lcm, prod
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from zlib import crc32

from flint import fmpz


def factors(n):
    result = [(int(p), int(e)) for p, e in fmpz(n).factor()]
    assert all(fmpz(p).is_prime() and e >= 1 for p, e in result)
    assert prod(p ** e for p, e in result) == n
    return result


def inputs(bound, minimum, maximum=None):
    atoms = factors(lcm(*range(1, bound + 1)))

    def visit(i, m, chosen):
        if maximum is not None and m > maximum:
            return
        if i == len(atoms):
            if m >= minimum:
                n = m + 1
                weight = [(p, p ** e) for p, e in factors(n)]
                data = [n, len(weight), *(v for pair in weight for v in pair),
                        len(chosen), *(v for pair in chosen for v in pair)]
                yield (" ".join(map(str, data)) + "\n").encode()
            return
        yield from visit(i + 1, m, chosen)
        p, exponent = atoms[i]
        for e in range(1, exponent + 1):
            w = p ** e
            yield from visit(i + 1, m * w, chosen + [(p, w)])

    yield from visit(0, 1, [])


def run(binary, bound, minimum, jobs, maximum=None, retained=False):
    digest = sha256()
    with TemporaryDirectory(dir=Path(binary).resolve().parent) as directory:
        paths = [Path(directory) / str(j) for j in range(jobs)]
        files = [path.open("wb") for path in paths]
        count = 0
        try:
            for line in inputs(bound, minimum, maximum):
                digest.update(line)
                files[crc32(line) % jobs].write(line)
                count += 1
        finally:
            for file in files:
                file.close()
        print(f"{count} weights, input SHA-256 {digest.hexdigest()}", flush=True)

        def search(path):
            with path.open("rb") as source:
                command = [str(Path(binary).resolve())] + (["--retained"] if retained else [])
                result = subprocess.run(command, stdin=source,
                                        capture_output=True, text=True, check=True)
            summary = re.fullmatch(r"(\d+) weights, (\d+) subsets, (\d+) smallest cuts, (\d+) rows\n",
                                   result.stderr)
            if not summary:
                raise RuntimeError(result.stderr)
            return tuple(map(int, summary.groups())), result.stdout.splitlines()

        with ThreadPoolExecutor(max_workers=jobs) as pool:
            outputs = list(pool.map(search, paths))
        totals = tuple(map(sum, zip(*(stats for stats, _ in outputs))))
        rows = sorted(line for _, batch in outputs for line in batch)
        assert totals[0] == count and totals[3] == len(rows)
        print("Weights, subsets, smallest cuts and full rows:", totals)
        for line in rows:
            print(line)
        return totals, rows


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("binary")
    parser.add_argument("--bound", type=int, default=53)
    parser.add_argument("--minimum-predecessor", type=int, default=10**9)
    parser.add_argument("--maximum-predecessor", type=int)
    parser.add_argument("--retained", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if not 2 <= args.bound <= 67 or args.minimum_predecessor < 2 or args.jobs < 1:
        parser.error("require 2 <= bound <= 67, minimum predecessor >= 2, and jobs >= 1")
    run(args.binary, args.bound, args.minimum_predecessor, args.jobs, args.maximum_predecessor, args.retained)
