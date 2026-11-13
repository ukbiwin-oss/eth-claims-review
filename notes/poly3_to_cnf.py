#!/usr/bin/env python3
"""Encode a sparse degree-3 polynomial over F2 as a 5-CNF whose satisfying assignments match the polynomial's
zeros one to one, and test that claim by brute force.

Status: unreviewed. Written 13 November 2026 to check one inference in section 4 of
`eth-claims-and-quantum-advantage.md` ("that polynomial encodes into a 5-CNF with about 700 variables whose
satisfying assignments match the polynomial's zeros one to one").

The polynomial
    f(x) = c + sum of monomials (mod 2),   each monomial a product of 1, 2 or 3 distinct variables of x_1..x_n,
is the kind Dalzell et al. (arXiv:1805.05224, section 2) attach to an IQP circuit: one Z, CZ or CCZ gate per
monomial of degree 1, 2 or 3, one qubit per variable. They take c = 0; c = 1 is supported here as well.

The encoding
  1. Split the monomials into groups G_1..G_C so that each group mentions at most 3 variables in all.
     (One monomial per group always works. Monomials on the same three variables can share a group.)
  2. Add chain variables s_1..s_{C-1}; s_j is the running parity of groups 1..j.
  3. Constraint j says   s_{j-1} + G_j(x) + s_j = 0 (mod 2),   with s_0 = c and s_C = 0 as constants.
     It mentions at most 2 chain variables and 3 x-variables: at most 5 variables.
  4. A constraint on v <= 5 variables becomes one clause of width v for each assignment of those variables that
     violates it. The result is a CNF with clauses of width at most 5.

Why the match is one to one: for a fixed x the constraints, read in order, force s_1, then s_2, ..., then
s_{C-1}; the last constraint then holds exactly when c + G_1(x) + ... + G_C(x) = f(x) = 0. So every zero of f
extends to exactly one satisfying assignment and every other x to none. The code does not rely on this
argument: the two checks below test it without using the chain structure.

Variable count: N = n + C - 1 <= n + M - 1 for M >= 1 monomials (N = n when there are none). Clause count: at
most 16 per constraint and 8 for the first and last, so at most 16M - 16 for M >= 2. These are figures for this
encoding, not lower bounds over all encodings.

Variables that occur in no monomial occur in no clause. They are free on both sides, so the match still holds,
but a model counter must take the variable count from the `p cnf N` header, not from the clauses.

Usage:
    python3 poly3_to_cnf.py test                 # run all checks and print the report (under a minute)
    python3 poly3_to_cnf.py count N M [SEED]     # variable and clause counts for random instances
    python3 poly3_to_cnf.py dimacs FILE          # FILE: first line "n c", then one monomial per line ("1 5 9")
"""
from __future__ import annotations

import itertools
import sys
import time

import numpy as np

Monomial = tuple[int, ...]          # sorted, distinct variable numbers, 1-based, length 1..3


# ----------------------------------------------------------------------------------------------------- encoding
def normalise(monomials) -> list[Monomial]:
    """Sort each monomial, reject bad ones, and cancel repeats (m + m = 0 over F2)."""
    seen: dict[Monomial, int] = {}
    for m in monomials:
        t = tuple(sorted(m))
        if not 1 <= len(t) <= 3 or len(set(t)) != len(t) or t[0] < 1:
            raise ValueError(f'not a monomial of degree 1-3 in distinct variables: {m!r}')
        seen[t] = seen.get(t, 0) ^ 1
    return [t for t, on in seen.items() if on]


def group(monomials: list[Monomial], greedy=True) -> list[list[Monomial]]:
    """Groups of monomials, each mentioning at most 3 variables. greedy=False: one monomial per group (the worst
    case, C = M). greedy=True: first fit, widest monomials first, so that a degree-3 monomial collects the
    lower-degree monomials inside its own three variables."""
    if not greedy:
        return [[m] for m in monomials]
    groups: list[tuple[set, list]] = []
    by_var: dict[int, list[int]] = {}
    for m in sorted(monomials, key=lambda t: (-len(t), t)):
        placed = False
        for gi in sorted({g for v in m for g in by_var.get(v, [])}):      # only groups that share a variable
            sup, members = groups[gi]
            if len(sup | set(m)) <= 3:
                members.append(m); new = set(m) - sup; sup |= new
                for v in new:
                    by_var.setdefault(v, []).append(gi)
                placed = True
                break
        if not placed:
            groups.append((set(m), [m]))
            for v in m:
                by_var.setdefault(v, []).append(len(groups) - 1)
    # second pass: merge any two groups whose variables together number at most 3 (for example two degree-1
    # monomials, or a degree-1 and a degree-2 monomial, on different variables)
    groups.sort(key=lambda g: -len(g[0]))
    merged: list[tuple[set, list]] = []
    for sup, members in groups:
        for msup, mmembers in merged:
            if len(msup | sup) <= 3:
                msup |= sup; mmembers += members
                break
        else:
            merged.append((set(sup), list(members)))
    return [members for _, members in merged]


def encode(n: int, monomials, const: int = 0, greedy=True):
    """-> (N, clauses, info). Variables 1..n are x; n+1..N are the chain variables. Clauses are lists of
    non-zero ints, DIMACS style (negative = negated)."""
    monos = normalise(monomials)
    if any(v > n for m in monos for v in m):
        raise ValueError('a monomial uses a variable above n')
    groups = group(monos, greedy)
    C = len(groups)
    if C == 0:                                   # f is the constant c: no clause if c = 0, the empty clause if c = 1
        return n, ([[]] if const & 1 else []), {'groups': 0, 'chain_variables': 0, 'monomials': 0, 'max_width': 0}
    N = n + C - 1
    clauses: list[list[int]] = []
    width = 0
    for j, members in enumerate(groups, 1):
        xs = sorted({v for m in members for v in m})
        prev = n + j - 1 if j > 1 else None      # s_{j-1}; for j = 1 it is the constant c
        nxt = n + j if j < C else None           # s_j;     for j = C it is the constant 0
        vs = xs + [v for v in (prev, nxt) if v is not None]
        width = max(width, len(vs))
        for bits in itertools.product((0, 1), repeat=len(vs)):
            val = dict(zip(vs, bits))
            parity = ((val[prev] if prev is not None else const & 1)
                      ^ (val[nxt] if nxt is not None else 0))
            for m in members:
                parity ^= all(val[v] for v in m)
            if parity:                           # this assignment violates the constraint: forbid it
                clauses.append([-v if b else v for v, b in zip(vs, bits)])
    return N, clauses, {'groups': C, 'chain_variables': C - 1, 'monomials': len(monos), 'max_width': width}


def dimacs(N: int, clauses) -> str:
    return f'p cnf {N} {len(clauses)}\n' + ''.join(' '.join(map(str, c)) + ' 0\n' for c in clauses)


# ------------------------------------------------------------------------------------------------- brute force
def poly_values(n: int, monomials, const: int = 0) -> np.ndarray:
    """f(x) for every x in {0,1}^n; x is the integer whose bit i-1 is x_i."""
    idx = np.arange(1 << n, dtype=np.uint32)
    f = np.full(1 << n, const & 1, dtype=np.uint8)
    for m in normalise(monomials):
        t = np.ones(1 << n, dtype=np.uint8)
        for v in m:
            t &= ((idx >> (v - 1)) & 1).astype(np.uint8)
        f ^= t
    return f


def measure(n: int, N: int, clauses, monomials) -> dict:
    """Bounds measured from the clause list itself (not the encoder's own report), and whether they hold."""
    M = len(normalise(monomials))
    width = max((len(c) for c in clauses), default=0)
    used = {abs(lit) for c in clauses for lit in c}
    ok = (width <= 5 and N <= n + max(0, M - 1)
          and len(clauses) <= (max(16 * M - 16, 8) if M else 1) and all(1 <= v <= N for v in used)
          and all(len({abs(l) for l in c}) == len(c) for c in clauses))
    return {'measured_width': width, 'bounds_ok': ok}


def check_exhaustive(n: int, monomials, const: int = 0, greedy=True, encoded=None, reference=None) -> dict:
    """Check 1 (no shortcut): evaluate the CNF on all 2^N assignments. Passes when, for every x, the number of
    satisfying assignments that extend x is 1 if f(x) = 0 and 0 if f(x) = 1. `encoded` substitutes another
    (N, clauses, info), used by the control tests; `reference` substitutes another table of f."""
    N, clauses, info = encoded or encode(n, monomials, const, greedy)
    if N > 24:
        raise ValueError(f'{N} CNF variables: too many for the exhaustive check')
    idx = np.arange(1 << N, dtype=np.uint32)
    bits = {v: ((idx >> (v - 1)) & 1).astype(bool) for v in range(1, N + 1)}
    sat = np.ones(1 << N, dtype=bool)
    for c in clauses:
        ok = np.zeros(1 << N, dtype=bool)
        for lit in c:
            ok |= bits[abs(lit)] if lit > 0 else ~bits[abs(lit)]
        sat &= ok
    del bits
    per_x = np.bincount(idx[sat] & ((1 << n) - 1), minlength=1 << n)
    f = poly_values(n, monomials, const) if reference is None else reference
    zeros = int((f == 0).sum())
    m = measure(n, N, clauses, monomials)
    return {'ok': bool(np.array_equal(per_x, (f == 0).astype(per_x.dtype))) and m['bounds_ok'], 'n': n, 'N': N,
            'clauses': len(clauses), 'zeros': zeros, 'satisfying': int(sat.sum()), **info, **m}


def check_propagation(n: int, monomials, const: int = 0, greedy=True, max_rounds=10_000, encoded=None,
                      sample: int = 0, rng=None) -> dict:
    """Check 2 (for more chain variables than check 1 can enumerate): for every x at once, fix x and run plain
    unit propagation on the clauses, in the order given, until nothing changes. Unit propagation only ever sets a
    variable to a value every extension must take. So if it ends with every variable set and no clause false,
    that is the one and only extension of x; if it ends with a clause false, x has no extension. Anything else
    (a variable left unset) counts as a failure of the check. The routine knows nothing about chains or parity:
    it sees only the clause list. With `sample` > 0 it runs from that many random x instead of all 2^n (for
    sizes where 2^n is out of reach) and f is evaluated directly on those x."""
    N, clauses, info = encoded or encode(n, monomials, const, greedy)
    if sample:
        size = sample
        bits = rng.integers(0, 2, size=(n + 1, size), dtype=np.int8)
    else:
        size = 1 << n
        idx = np.arange(size, dtype=np.uint32)
        bits = np.zeros((n + 1, size), dtype=np.int8)
        for v in range(1, n + 1):
            bits[v] = (idx >> (v - 1)) & 1
    val = np.full((N + 1, size), -1, dtype=np.int8)               # -1 unset, 0 false, 1 true
    val[1:n + 1] = bits[1:]
    conflict = np.zeros(size, dtype=bool)
    for rounds in range(1, max_rounds + 1):
        changed = False
        for c in clauses:
            satisfied = np.zeros(size, dtype=bool)
            unset = np.zeros(size, dtype=np.int8)
            for lit in c:
                col = val[abs(lit)]
                satisfied |= col == (1 if lit > 0 else 0)
                unset += col == -1
            conflict |= ~satisfied & (unset == 0)
            unit = ~satisfied & (unset == 1)
            if unit.any():
                for lit in c:
                    col = val[abs(lit)]
                    put = unit & (col == -1)
                    if put.any():
                        col[put] = 1 if lit > 0 else 0
                        changed = True
        if not changed:
            break
    complete = (val[1:] != -1).all(axis=0)
    if sample:
        f = np.full(size, const & 1, dtype=np.uint8)
        for m in monomials:                                        # the raw list: repeats cancel by themselves
            t = np.ones(size, dtype=np.uint8)
            for v in m:
                t &= bits[v].astype(np.uint8)
            f ^= t
    else:
        f = poly_values(n, monomials, const)
    one = complete & ~conflict                                     # exactly one extension
    none = conflict                                                # no extension
    undecided = int((~one & ~none).sum())
    m = measure(n, N, clauses, monomials)
    ok = (undecided == 0 and bool(np.array_equal(one, f == 0)) and bool(np.array_equal(none, f == 1))
          and m['bounds_ok'])
    return {'ok': ok, 'n': n, 'N': N, 'clauses': len(clauses), 'zeros': int((f == 0).sum()),
            'unique_extension': int(one.sum()), 'no_extension': int(none.sum()), 'undecided': undecided,
            'rounds': rounds, **info, **m}


def naive_values(n: int, raw_monomials, const: int = 0) -> np.ndarray:
    """f(x) for every x, by a second, plain-Python route that shares no code with `normalise` or `poly_values`:
    each listed monomial is evaluated and added mod 2, repeats included."""
    out = np.zeros(1 << n, dtype=np.uint8)
    for x in range(1 << n):
        v = const & 1
        for m in raw_monomials:
            v ^= all(x >> (i - 1) & 1 for i in m)
        out[x] = v
    return out


# ------------------------------------------------------------------------------------------------- instances
def random_poly(rng, n: int, M: int, degrees=(1, 2, 3), weights=None) -> list[Monomial]:
    """M distinct random monomials; each one's degree is drawn from `degrees` (with `weights`)."""
    available = sum(len(list(itertools.combinations(range(n), d))) for d in set(degrees)) if n <= 12 else M
    if M > available:
        raise ValueError(f'only {available} distinct monomials of these degrees exist on {n} variables')
    out: set[Monomial] = set()
    while len(out) < M:
        d = int(rng.choice(degrees, p=weights))
        out.add(tuple(sorted(int(v) + 1 for v in rng.choice(n, size=d, replace=False))))
    return sorted(out)


def all_monomials(n: int) -> list[Monomial]:
    return [m for d in (1, 2, 3) for m in itertools.combinations(range(1, n + 1), d)]


# ------------------------------------------------------------------------------------------------------- report
def run_tests() -> int:
    rng = np.random.default_rng(20261113)
    t0 = time.time()
    failures = 0
    print('poly3_to_cnf.py test report')
    print(f'python {sys.version.split()[0]}, numpy {np.__version__}, seed 20261113')
    print('Every PASS also requires the bounds: width <= 5 and clauses <= 16M - 16 for M >= 2 (8 for M = 1), measured '
          'from the clause list; N <= n + M - 1 for M >= 1.\n')

    # A. every polynomial on 3 variables (7 monomials -> 128 polynomials), both constant terms, both groupings
    cases = bad = 0
    monos3 = all_monomials(3)
    for mask in range(1 << len(monos3)):
        poly = [m for i, m in enumerate(monos3) if mask >> i & 1]
        for const in (0, 1):
            for greedy in (True, False):
                r = check_exhaustive(3, poly, const, greedy); cases += 1; bad += not r['ok']
    print(f'A. all 128 polynomials of degree at most 3 on 3 variables x constant term 0/1 x both groupings (greedy gives one group '
          f'here, so only the single grouping builds chains), every CNF assignment enumerated: '
          f'{cases - bad}/{cases} pass')
    failures += bad

    # B. every polynomial on 4 variables (14 monomials -> 16,384), both constant terms, greedy grouping
    cases = bad = 0
    monos4 = all_monomials(4)
    worstN = 0
    for mask in range(1 << len(monos4)):
        poly = [m for i, m in enumerate(monos4) if mask >> i & 1]
        for const in (0, 1):
            r = check_exhaustive(4, poly, const, True); cases += 1; bad += not r['ok']; worstN = max(worstN, r['N'])
    print(f'B. all 16,384 polynomials of degree at most 3 on 4 variables x constant term 0/1, greedy grouping (chains of up to '
          f'{worstN - 4} chain variables), every CNF assignment enumerated: {cases - bad}/{cases} pass')
    failures += bad

    # C. random sparse polynomials, exhaustive
    print('\nC. random polynomials, every CNF assignment enumerated (check 1)')
    print('   n   M  c  grouping   N  clauses  width      zeros  satisfying  result')
    for n, M, greedy in [(5, 8, False), (6, 10, False), (8, 12, False), (10, 12, False), (12, 12, False),
                         (14, 10, False), (16, 9, False), (18, 7, False), (20, 5, False),
                         (6, 19, True), (8, 17, True), (10, 15, True), (12, 13, True), (16, 9, True)]:
        for rep in range(3):
            poly = random_poly(rng, n, M); const = rep % 2 if rep < 2 else int(rng.integers(2))
            r = check_exhaustive(n, poly, const, greedy)
            failures += not r['ok']
            print(f"  {n:2d} {M:3d}  {const}  {'greedy' if greedy else 'single':8s} {r['N']:3d} {r['clauses']:8d} "
                  f"{r['measured_width']:6d} {r['zeros']:10d} {r['satisfying']:11d}  {'PASS' if r['ok'] else 'FAIL'}")

    # D. longer chains at up to 20 x-variables: unit propagation from every x (check 2)
    print('\nD. random polynomials, unit propagation from every x (check 2)')
    print('   n   M  c  grouping   N  clauses  width      zeros  unique-ext  no-ext  undecided  result')
    for n, M, greedy in [(10, 30, False), (12, 40, False), (14, 40, False), (16, 40, False), (18, 40, False),
                         (20, 40, False), (20, 60, False), (20, 60, True), (16, 100, True), (16, 100, False)]:
        poly = random_poly(rng, n, M); const = int(rng.integers(2))
        r = check_propagation(n, poly, const, greedy)
        failures += not r['ok']
        print(f"  {n:2d} {M:3d}  {const}  {'greedy' if greedy else 'single':8s} {r['N']:3d} {r['clauses']:8d} "
              f"{r['measured_width']:6d} {r['zeros']:10d} {r['unique_extension']:11d} {r['no_extension']:7d} "
              f"{r['undecided']:10d}  {'PASS' if r['ok'] else 'FAIL'}")

    # E. a second, separately written evaluation of f, and repeated monomials
    cases = bad = 0
    for _ in range(300):
        n = int(rng.integers(3, 9)); M = int(rng.integers(1, 10))
        poly = random_poly(rng, n, min(M, len(all_monomials(n))))
        raw = poly + [poly[int(i)] for i in rng.integers(0, len(poly), size=int(rng.integers(0, 4)))]   # repeats
        const = int(rng.integers(2)); greedy = bool(rng.integers(2))
        r = check_exhaustive(n, raw, const, greedy, reference=naive_values(n, raw, const))
        cases += 1; bad += not r['ok']
    print(f'\nE. 300 random polynomials with up to 3 repeated monomials each, zeros taken from a second plain-Python evaluation '
          f'of the raw list (shares no code with the encoder): {cases - bad}/{cases} pass')
    failures += bad

    # F. controls: the checks themselves must reject broken encodings
    # A mutant can still be a correct encoding (the altered clause concerned a combination of values that cannot
    # occur), so check 1, which enumerates everything, is the judge of that. Required: every mutant check 1
    # rejects is also rejected by check 2, and a wrong constant term is rejected by both every time.
    n_mut = wrong1 = wrong2 = missed = const_total = const_caught = 0
    for trial in range(90):
        n = int(rng.integers(5, 9)); poly = random_poly(rng, n, int(rng.integers(4, 9))); const = int(rng.integers(2))
        N, clauses, info = encode(n, poly, const, False)
        k = int(rng.integers(len(clauses)))
        kind = trial % 3
        if kind == 0:                                              # drop a clause
            broken = clauses[:k] + clauses[k + 1:]
        elif kind == 1:                                            # flip one literal
            c = list(clauses[k]); c[0] = -c[0]; broken = clauses[:k] + [c] + clauses[k + 1:]
        else:                                                      # wrong constant term
            broken = encode(n, poly, const ^ 1, False)[1]
        a = check_exhaustive(n, poly, const, False, encoded=(N, broken, info))
        b = check_propagation(n, poly, const, False, encoded=(N, broken, info))
        n_mut += 1; wrong1 += not a['ok']; wrong2 += not b['ok']; missed += (not a['ok']) and b['ok']
        if kind == 2:
            const_total += 1; const_caught += (not a['ok']) and (not b['ok'])
    print(f'F. controls: {n_mut} altered encodings (a clause dropped, a literal flipped, or the wrong constant term). '
          f'Check 1 rejects {wrong1}; the other {n_mut - wrong1} are still correct encodings. Check 2 rejects {wrong2}, '
          f'including every one check 1 rejects ({missed} missed). Wrong constant term: {const_caught}/{const_total} '
          f'rejected by both.')
    failures += missed + (const_total - const_caught)

    # G. the target size: 208 variables, 500 monomials
    print('\nG. 208 variables and 500 monomials: unit propagation from 5,000 random x (2^208 cannot be enumerated)')
    print('   mix                           c  grouping   N  clauses  width  f=0  unique-ext  f=1  no-ext  undecided  result')
    for name, degrees in [('all degree 3 (CCZ only)', (3,)), ('degrees 1, 2, 3 (Z, CZ, CCZ)', (1, 2, 3))]:
        for greedy in (False, True):
            poly = random_poly(rng, 208, 500, degrees); const = int(rng.integers(2))
            r = check_propagation(208, poly, const, greedy, sample=5000, rng=rng)
            failures += not r['ok']
            print(f"   {name:29s} {const}  {'greedy' if greedy else 'single':8s} {r['N']:3d} {r['clauses']:7d} "
                  f"{r['measured_width']:6d} {r['zeros']:5d} {r['unique_extension']:10d} "
                  f"{5000 - r['zeros']:5d} {r['no_extension']:6d} {r['undecided']:9d}  {'PASS' if r['ok'] else 'FAIL'}")

    print('\nVariable counts for 208 variables and 500 monomials (the 500-gate IQP case of the report)')
    count_report(208, 500, rng)
    print(f'\n{"ALL CHECKS PASS" if not failures else f"{failures} FAILURE(S)"}  ({time.time() - t0:.0f} s)')
    return 1 if failures else 0


def count_report(n: int, M: int, rng, reps=20) -> None:
    print(f'  worst case for this encoding, any polynomial: n + M - 1 = {n + M - 1} variables (one chain variable per '
          f'monomial but one); at most 16M - 16 = {16 * M - 16} clauses; width at most 5')
    print('  this greedy grouping on uniformly random monomials (20 draws each; min / median / max):')
    mixes = [('all degree 3 (CCZ only)', (3,), None),
             ('degrees 1, 2, 3 equally likely (Z, CZ, CCZ)', (1, 2, 3), None),
             ('degrees 2 and 3 equally likely (CZ, CCZ)', (2, 3), None),
             ('degree 2 only (CZ only; quadratic, so easy anyway)', (2,), None)]
    for name, degrees, w in mixes:
        Ns, Cs, ws = [], [], []
        for _ in range(reps):
            N, clauses, info = encode(n, random_poly(rng, n, M, degrees, w), 0, True)
            Ns.append(N); Cs.append(len(clauses)); ws.append(info['max_width'])
        Ns.sort(); Cs.sort()
        print(f'    {name:52s} variables {Ns[0]} / {Ns[len(Ns) // 2]} / {Ns[-1]}   clauses {Cs[0]} / '
              f'{Cs[len(Cs) // 2]} / {Cs[-1]}   width {max(ws)}')


def main(argv) -> int:
    if len(argv) >= 2 and argv[1] == 'test':
        return run_tests()
    if len(argv) >= 4 and argv[1] == 'count':
        count_report(int(argv[2]), int(argv[3]), np.random.default_rng(int(argv[4]) if len(argv) > 4 else 0))
        return 0
    if len(argv) == 3 and argv[1] == 'dimacs':
        lines = [ln.split() for ln in open(argv[2]) if ln.strip() and not ln.startswith('#')]
        n, const = int(lines[0][0]), int(lines[0][1]) if len(lines[0]) > 1 else 0
        N, clauses, info = encode(n, [tuple(map(int, ln)) for ln in lines[1:]], const)
        sys.stdout.write(f'c {info}\n' + dimacs(N, clauses))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
