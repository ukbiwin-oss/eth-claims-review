#!/usr/bin/env python3
"""Brute-force check of steps 1 and 2 of conjecture7.md (unreviewed; 13 November 2026).

For a clause list f on n variables with m clauses of at most 3 literals:
  step 1   gap(f) = sum over x of (-1)^f(x) = 2^n - 2S,          S  = number of satisfying assignments of f
  step 2   n' <= 3m  and  gap(f) = 2^(n-n') * gap(f'),            f' = f restricted to the n' variables it uses

gap(f) is computed directly from its definition over all 2^n assignments; the right-hand sides are computed
from a separate count over the 2^n' assignments of the used variables only. This tests the arithmetic, not the
running-time claim.

Usage: python3 gap_check.py
"""
import itertools
import random
import sys


def evaluate(clauses, assignment) -> int:
    """1 if every clause has a true literal. A literal is (variable index, wanted value)."""
    return int(all(any(assignment[v] == want for v, want in c) for c in clauses))


def gap_direct(n, clauses) -> int:
    return sum((-1) ** evaluate(clauses, x) for x in itertools.product((0, 1), repeat=n))


def reduce(clauses):
    """Drop unused variables and rename the used ones 0..n'-1."""
    used = sorted({v for c in clauses for v, _ in c})
    name = {v: i for i, v in enumerate(used)}
    return len(used), [[(name[v], want) for v, want in c] for c in clauses]


def count(n, clauses) -> int:
    return sum(evaluate(clauses, x) for x in itertools.product((0, 1), repeat=n))


def check(n, clauses) -> bool:
    m = len(clauses)
    g = gap_direct(n, clauses)
    step1 = g == 2 ** n - 2 * count(n, clauses)
    n2, reduced = reduce(clauses)
    s2 = count(n2, reduced)
    step2 = n2 <= 3 * m and g == 2 ** (n - n2) * (2 ** n2 - 2 * s2)
    same_answer = (g != 0) == (2 * s2 != 2 ** n2)
    return step1 and step2 and same_answer


def main() -> int:
    rng = random.Random(20261113)
    cases = bad = balanced = with_unused = 0
    edge = {'no clauses, n=0': (0, []), 'no clauses, n=5': (5, []), 'one empty clause, n=4': (4, [[]]),
            'x1 alone, n=1 (balanced)': (1, [[(0, 1)]]), 'x1 alone, n=6 (balanced)': (6, [[(0, 1)]]),
            'x1 and not x1, n=3': (3, [[(0, 1)], [(0, 0)]])}
    for name, (n, clauses) in edge.items():
        ok = check(n, clauses); cases += 1; bad += not ok
        print(f'edge case: {name:28s} gap = {gap_direct(n, clauses):4d}  {"PASS" if ok else "FAIL"}')
    for _ in range(20000):
        n = rng.randint(1, 12)
        m = rng.randint(0, 8)
        pool = rng.sample(range(n), rng.randint(1, n))            # the clauses draw on a subset: some variables unused
        clauses = []
        for _ in range(m):
            width = rng.randint(1, min(3, len(pool)))
            clauses.append([(v, rng.randint(0, 1)) for v in rng.sample(pool, width)])
        ok = check(n, clauses); cases += 1; bad += not ok
        balanced += gap_direct(n, clauses) == 0
        with_unused += reduce(clauses)[0] < n
    print(f'\n{cases - bad}/{cases} cases pass (20,000 random clause lists, n up to 12, m up to 8, clause width 1-3; '
          f'{with_unused} had unused variables, {balanced} had gap = 0)')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
