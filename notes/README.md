# Notes: the physics report's inferences, written out and tested

These notes take inferences that [the physics report](../eth-claims-and-quantum-advantage.md) marks as unreviewed and make each one small enough to check. **They are still unreviewed.** No expert has read them. Each note says what was checked and what was not. Everything here is conditional on the October 2026 claims being correct; the notes take no view on that.

| Note | Inference in the report | Result |
|---|---|---|
| [conjecture7.md](conjecture7.md) | An exact #3-SAT counter in 2^{o(n)} time refutes Morimae–Tamaki Conjecture 7 (section 4, rows A7 and A8) | Written out in five steps with its assumptions stated. No error found so far. One addition: Theorem 8 of the paper rests on the same conjecture as Theorem 7. |
| [poly3_to_cnf.py](poly3_to_cnf.py) | A 500-gate IQP circuit on 208 qubits encodes into a 5-CNF with "about 700" variables whose satisfying assignments match the polynomial's zeros one to one (section 4, row A12b) | The encoding is proved in the file's header and tested. For this encoding the worst case is 707 variables and 7,984 clauses of width at most 5. |

## Files

- `conjecture7.md`: the argument, the points an expert should check, and what was read.
- `gap_check.py`, `gap_check_output.txt`: a sanity check of the note's arithmetic on 20,006 small clause lists.
- `poly3_to_cnf.py`: the encoder and its tests. `python3 poly3_to_cnf.py test` reproduces the output file in under a minute.
- `poly3_to_cnf_test_output.txt`: the test run of 13 November 2026.

## The encoder

The correctness of the encoding rests on the short argument in the header of `poly3_to_cnf.py`: for a fixed x the constraints force each chain variable in turn, and the last constraint holds exactly when f(x) = 0. The tests are there to catch a mistake in the code, and they stop well short of the target size.

Two checks are used. They share the encoder but not the argument: neither knows about chains or parity.

1. **Every assignment enumerated.** For CNFs of up to 24 variables, all 2^N assignments are evaluated. For every x the number of satisfying assignments above it must be 1 when f(x) = 0 and 0 when f(x) = 1.
2. **Unit propagation.** x is fixed and plain unit propagation is run on the clause list. It must end with every variable forced and no clause false when f(x) = 0, and with a false clause when f(x) = 1. This is a sufficient test: it can fail on a correct CNF, but when it passes the extension is unique.

What was run (all pass):

- Check 1 on every polynomial in 3 variables and every polynomial in 4 variables, with both constant terms, and on 42 random polynomials with 5 to 20 variables.
- Check 2 from every x on 10 random polynomials with 10 to 20 variables and up to 100 monomials. The longest chain has 99 chain variables (16 variables, 100 monomials).
- Check 1 against a second, separately written evaluation of the polynomial, on 300 random polynomials that include repeated monomials.
- Controls: 90 deliberately altered encodings. Check 1 rejects 75; the other 15 are still correct encodings, because the altered clause concerned values that cannot occur. Check 2 rejects every one that check 1 rejects, and a wrong constant term is rejected by both every time.
- At the target size, 208 variables and 500 monomials, check 2 from 5,000 random x on four instances. 2^208 assignments cannot be enumerated, so this is a sample, not a verification.
- Every pass also requires the bounds measured from the clause list itself: width at most 5, at most n + M − 1 variables, at most 16M − 16 clauses.

## The count for 208 variables and 500 monomials

- **707 variables and 7,984 clauses** when every monomial has degree 3 on its own triple of variables (500 CCZ gates). This is the worst case for this encoding, n + M − 1 and 16M − 16. It is not a lower bound over all encodings.
- Fewer when lower-degree monomials can share a constraint. With this file's grouping heuristic on uniformly random gates, a mix of Z, CZ and CCZ gave 488 to 518 variables over 20 draws, and CZ with CCZ gave 577 to 609.

So the report's "about 700" is right for the case it describes.

## What this does and does not establish

It establishes the encoding and the count. It does not establish anything about running time, which is the part of the report's inference that matters: that a counter running in 2^{o(N)} time on such CNFs would simulate circuits with O(n) gates in 2^{o(n)} time, where the 500-gate figures assume a cost proportional to the number of gates times 2^{an}. That step is unchanged and is still an inference.

Two limits are worth stating. At this size N is about 3.4n, so a counter running in 2^{δN} time beats plain enumeration of the 2^208 inputs only if δ is below 208/707, about 0.29. And the encoder says nothing useful about dense instances, where the number of monomials grows like n³ and this encoding gives far more CNF variables than the polynomial has.

## How these notes were checked

By the tests above, by reading the cited passages of the two source papers, and by asking a language model to referee a draft of each note. The model found no error in the five steps of note 1 or in the encoding; it did find a wrong justification in note 1 and several weak tests and loose statements here, which were corrected. That is proofreading. It is not review by an expert, and the notes stay marked unreviewed until one has read them.

## Not done yet

- The QAOA case, which the report says "was not worked through".
- How fast the clause width k may grow before each lab's exponent stops being o(n) (row A2).
- A Lean statement of note 1.
