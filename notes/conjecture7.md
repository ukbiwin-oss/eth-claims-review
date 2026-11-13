# Note 1. An exact #3-SAT counter in 2^{o(n)} time refutes Morimae–Tamaki Conjecture 7

**Status: unreviewed inference.** Written 13 November 2026. No expert has checked it, and we found no published source that states it. It is conditional throughout: it says what follows *if* the October 2026 counting claims are correct, and takes no view on whether they are. It expands the bracketed inference in section 4 of [the physics report](../eth-claims-and-quantum-advantage.md) (table rows A7 and A8). The implication is short; the note exists so that each step can be accepted or rejected quickly.

## The conjecture

Source: T. Morimae and S. Tamaki, "Fine-grained quantum computational supremacy", [arXiv:1901.01637](https://arxiv.org/abs/1901.01637) (version 4 read; Quantum Information and Computation 19, 1089–1115).

For a Boolean function f : {0,1}^n → {0,1}, equation (7) of the paper defines

    gap(f) = Σ_{x ∈ {0,1}^n} (−1)^{f(x)}.

**Conjecture 7** (paraphrased closely). Any non-deterministic algorithm that decides whether gap(f) ≠ 0 or gap(f) = 0, given a description of a 3-CNF f with m clauses, "needs 2^{Ω(m)} time".

The bound is in the number of clauses m, not the number of variables n. The paper's proof of its Theorem 7 (Appendix H) uses the conjecture in exactly this form: from a 3-CNF with m clauses it builds a circuit whose acceptance probability is non-zero exactly when gap(f) ≠ 0, derives a decision "in non-deterministic 2^{o(m)} time", and calls that a contradiction with Conjecture 7. Time is counted in m alone there.

## The hypothesis

**(H)** There are one deterministic algorithm, one constant a ≥ 0 and one function g with g(n) finite for every n ≥ 0 and g(n) = o(n), such that for every 3-CNF φ with n variables and m clauses the algorithm outputs the exact number #φ of satisfying assignments within

    (n + m + 1)^a · 2^{g(n)}   steps, in the worst case.

"3-CNF" here means a list of clauses with at most 3 literals each, degenerate cases included (repeated clauses, a clause containing a variable and its negation, an empty clause, no clauses).

Both labs claim a deterministic, worst-case, exact counter with g(n) = O(n log log n / √log n); see the comparison table in [the P vs NP report](../p-vs-np-and-eth-claims.md). Only the properties listed in (H) are used below.

## Convention on input length

**(E)** The input has length at most 2^{o(m)}; for example, variable names are numbers of at most poly(m) bits.

The conjecture counts time in m, so some such convention is needed for it to make sense, and the paper's Appendix H makes the same one implicitly when it counts the cost of building a circuit from f in m alone. The algorithm below never needs the value of n.

## The argument

Let f be a 3-CNF on variables x_1..x_n with m clauses. Write S for the number of x with f(x) = 1.

**Step 1. The gap is a fixed function of the count.** The sum has S terms equal to −1 and 2^n − S terms equal to +1, so

    gap(f) = 2^n − 2S,   and   gap(f) = 0  ⇔  S = 2^{n−1}.

**Step 2. Unused variables can be dropped, so n ≤ 3m.** Let n′ be the number of variables that occur in at least one clause. Each clause names at most 3 variables, so n′ ≤ 3m. Let f′ be the same list of clauses read as a function of those n′ variables, with S′ satisfying assignments. A variable in no clause is free, so every satisfying assignment of f′ extends in 2^{n−n′} ways:

    S = 2^{n−n′} · S′,   hence   gap(f) = 2^n − 2S = 2^{n−n′} · (2^{n′} − 2S′) = 2^{n−n′} · gap(f′).

The factor 2^{n−n′} is not zero, so gap(f) ≠ 0 ⇔ gap(f′) ≠ 0. Finding the used variables and renaming them 1..n′ takes time polynomial in the input length. (Edge cases: with no clauses, f′ is the constant 1 on the single empty assignment and gap(f′) = −1; an empty clause makes S′ = 0 and gap(f′) = 2^{n′}. Both are non-zero and both follow from the same formula.)

**Step 3. Run the counter on f′.** By (H) this returns S′ exactly within (n′ + m + 1)^a · 2^{g(n′)} steps. Put G(N) = max over k ≤ N of g(k), which is finite because g is. Then G is non-decreasing and still o(N): given ε > 0, g(k) ≤ εk for all k beyond some K, so G(N) ≤ max(G(K), εN), which is at most εN once N is large. With n′ ≤ 3m,

    steps ≤ (4m + 1)^a · 2^{G(3m)} = 2^{a·log₂(4m + 1) + G(3m)} = 2^{o(m)}.

The step through G matters: n′ may be much smaller than 3m, and g need not be monotone.

**Step 4. Compare.** Compute 2S′ and 2^{n′}, two integers of at most n′ + 2 bits, and answer "gap ≠ 0" exactly when they differ. This costs time polynomial in m. By steps 1 and 2 the answer is correct for f.

**Step 5. Deterministic is a special case of non-deterministic.** A deterministic algorithm is a non-deterministic algorithm with one computation path, which accepts exactly when that path accepts. A non-deterministic decision procedure for this problem is one-sided, and Appendix H fixes the side: accept exactly when gap(f) ≠ 0, as in step 4. Deterministic time is closed under complement, so the other orientation is covered too. Steps 2 to 4 are therefore a non-deterministic algorithm that decides gap(f) ≠ 0 against gap(f) = 0 for every 3-CNF with m clauses.

**Conclusion.** Under (E), the total time is polynomial in the input length plus 2^{o(m)}, which is 2^{o(m)}. For every constant c > 0 that is below 2^{cm} on all inputs with m large enough. This contradicts "needs 2^{Ω(m)} time" under each reading of the conjecture: a constant for each algorithm or one constant for all, and a bound holding for all large m or only for infinitely many. So:

> If (H) holds, then (under the input convention E) Morimae–Tamaki Conjecture 7 is false.

Put the other way round, Conjecture 7 implies a counting form of the exponential-time hypothesis stated in m: no single deterministic algorithm counts 3-CNF solutions in 2^{o(m)} time. Hypotheses of this kind go under the name #ETH (Dell et al., [arXiv:1206.1775](https://arxiv.org/abs/1206.1775)).

## What follows, and what does not

- **Theorems 7 and 8 of the paper lose their premise.** Both begin "Assume that Conjecture 7 is true": Theorem 7 for Clifford+T circuits with t T gates, Theorem 8 for t-qubit Pauli-based computation, each a 2^{o(t)} lower bound on sampling within multiplicative error. They are not shown false, and (H) does not give a 2^{o(t)} classical sampler for those circuits; the theorems become statements with a false hypothesis. The physics report lists Theorem 7 only (row A8); Theorem 8 has the same premise and belongs in the same row.
- **This argument needs more than a decision algorithm.** It uses a test of whether S′ = 2^{n′−1}, which the exact count gives. A 2^{o(n)} algorithm that only decides satisfiability is not known to give that test. Whether a refutation of ETH alone would reach Conjecture 7 by another route is open; this note does not exclude it.
- **A family of algorithms is not enough either.** (H) asks for one algorithm with one sub-linear exponent. A separate algorithm for each δ > 0 running in time 2^{δn}, which is the usual form of "#ETH is false", would not contradict the conjecture as literally stated or as used in Appendix H.
- **No sparsification is needed.** Arguments that pass from n to m for ETH use the sparsification lemma. Here the conjecture is already stated in m and n′ ≤ 3m is immediate.
- **Conjectures 4 to 6 of the paper are outside this note.** They concern general, logarithmic-depth and reversible circuits, with bounds of the form 2^{(1−a)n}. Nothing here addresses them; the physics report discusses them separately.
- **The paper's own motivation.** It says the only known way to decide gap(f) ≠ 0 is to solve #SAT, and that Conjecture 7 "could be suggested by ETH". That is a heuristic remark, not an implication. What this note shows is the reverse direction only: a fast enough exact counter removes the conjecture.

## Points an expert should check

1. **Reading of "needs 2^{Ω(m)} time".** Taken to mean: no algorithm of the stated kind runs in 2^{o(m)} time in the worst case. The paper's Appendix H uses it this way.
2. **Machine model.** Conjecture 7 says "non-deterministic algorithm" with no model named. The labs state their bounds in bit operations and on Turing machines. Standard models simulate one another with polynomial overhead, which 2^{o(m)} absorbs. If the conjecture were meant for a model in which that fails, step 5 would need restating.
3. **Input length.** The conclusion depends on convention (E). Without some bound on the input length in terms of m, neither the conjecture nor its refutation is well posed.
4. **Clause width and degenerate clauses.** (H) is stated for clauses of at most 3 literals, degenerate cases included. If a lab's counter is stated only for clauses of exactly 3 distinct variables, pad first: with two fresh variables z₁, z₂ shared by all clauses, (a ∨ b) becomes (a ∨ b ∨ z₁)(a ∨ b ∨ ¬z₁) and a single literal is padded twice. The count is multiplied by a known power of 2 (at most 4), m by at most 4, and n′ grows by at most 2, so the argument is unchanged; empty clauses and m = 0 are answered directly.
5. **Uniformity in (H).** (H) needs one algorithm, one constant a and one function g for all 3-CNFs. Whether each lab's formal statement gives exactly that is a question about the statements, discussed in the P vs NP report; this note does not settle it.
6. **Determinism.** A randomised counter with bounded error would not do for step 5; a zero-error one would. The labs' claims are taken to be deterministic, as the reports record.

## What was checked and what was not

Checked (13 November 2026): the wording of Conjecture 7, equation (7), Theorems 7 and 8 and the proof of Theorem 7, read from the text of arXiv:1901.01637 version 4. A language model was also asked to referee a draft of this note; it found no error in steps 1 to 5 and asked for the conventions now stated in (H), (E) and step 5. That is a proofreading pass, not a review by an expert.

`gap_check.py` runs the identities of steps 1 and 2 on 20,006 small clause lists as a sanity check on the arithmetic and the edge cases. The identities have two-line proofs, given above; the script adds no weight beyond catching slips.

Not checked: the journal version; whether the conjecture numbering differs in arXiv versions 1 to 3; the labs' theorem statements, which are taken from the two reports in this repository, not re-read; the Dell et al. paper beyond its role as the usual reference for #ETH; whether (H) is true.
