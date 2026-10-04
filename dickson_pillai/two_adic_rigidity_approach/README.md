# 2-Adic Defect Rigidity, Carry Automata, and Dynamical Repulsion

## 1. Overview & Document Metadata
* **Title:** *On the 2-Adic Rigidity, Carry Automata, and Dynamical Repulsion of the Dickson–Pillai Defect in Waring's Problem*
* **Target Publication:** *Journal of Number Theory* / *Acta Arithmetica* / *Mathematics of Computation*
* **Primary MSC (2020):** 11P05 (Waring's problem and variants), 11J86 (Linear forms in logarithms), 11S80 ($p$-adic analytic methods), 37B10 (Symbolic dynamics)
* **Codebase & Formal Verification:** Lean 4 formalization (Mathlib-compatible, zero non-standard axioms) in `formalization/TwoAdicDefectRigidity.lean`.

---

## 2. Core Proven Structural Theorems

1. **Theorem 1 (Residue Exclusion Modulo 8):**
   For all $k \ge 3$:
   $$D_k \bmod 8 = \begin{cases} 7 & \text{if } k \text{ is even}, \\ 5 & \text{if } k \text{ is odd}. \end{cases}$$
   * Immediate Corollaries: $D_k$ is strictly odd; $D_k \notin \{1, 3\}$ for all $k \ge 3$; exactly $75\%$ of all residue classes modulo $8$ are unconditionally eliminated.

2. **Theorem 2 (Exact 2-Adic Valuation for Even Indices):**
   By Lifting the Exponent Lemma (LTE) on $9^m - 1 = 8 V_m$, for all even $k \ge 6$:
   $$v_2(D_k + 1) = v_2(k) + 2 \iff D_k = 2^{v_2(k)+2} W - 1 \quad (W \text{ odd}).$$

3. **Theorem 3 (Exact 2-Adic Valuation for Odd Indices):**
   * For all $k \equiv 3 \pmod 4$ ($k \ge 5$): $D_k \equiv 5 \pmod{16} \iff v_2(D_k + 3) = 3$.
   * For all $k \equiv 1 \pmod 4$ ($k - v_2(k-1) \ge 3$): $v_2(D_k + 3) = v_2(k-1) + 2$.

4. **Theorem 4 (Exact 4-State Carry Drift Recurrence):**
   For all $k \ge 1$:
   $$3 D_k - D_{k+1} = C_k 2^k \iff D_{k+1} = 3 D_k - C_k 2^k,$$
   where $C_k = 3 q_k - 2 q_{k+1} + 1 \in \{-1, 0, 1, 2\}$ is governed by the cylinder partition on $(q_k \bmod 2, \delta_k)$:
   * $q_k$ even, $\delta_k \in [0, 2/3) \implies C_k = 1$
   * $q_k$ even, $\delta_k \in [2/3, 1) \implies C_k = -1$
   * $q_k$ odd, $\delta_k \in [0, 1/3) \implies C_k = 2$
   * $q_k$ odd, $\delta_k \in [1/3, 1) \implies C_k = 0$

5. **Theorem 5 (State Collapse Under Failure):**
   Any hypothetical failure $D_k < q_k$ forces:
   $$\delta_k > 1 - \left(\frac{3}{4}\right)^k > \frac{2}{3} \quad (\forall k \ge 5).$$
   Consequently, states $C_k \in \{1, 2\}$ are **strictly forbidden**, forcing $C_k \in \{-1, 0\}$.

6. **Theorem 6 (Dynamical Repulsion & Strict Failure Isolation):**
   Failures cannot cluster: if $k \ge 5$ is a failure index, step $k+1$ unconditionally satisfies:
   $$D_{k+1} > 2 q_{k+1}.$$
   Any failure orbit exponentially escapes via Branch 2 doubling or Branch 1 overshoot above $2^k$.

7. **Theorem 7 (Universal 2-Adic Fixed-Point Reduction):**
   Every even failure index $k = 2^{v_2(k)} u \ge 6$ satisfies an exact non-linear 2-adic fixed-point identity modulo $2^{k - v_2(k) - 2}$ in terms of the universal constant $\nu := \frac{\ln_2 3}{4} \in \mathbb{Z}_2^\times$:
   $$u \nu + W + \sum_{j=2}^\infty \frac{2^{(j-1)(v_2(k)+2)}}{j} W^j \equiv 0 \pmod{2^{k - v_2(k) - 2}}.$$
   The first-order linear approximation holds modulo $2^{v_2(k)+1}$, which allows linear Wang--Guy--Davenport reconstruction when $v_2(k) \ge \frac{1}{2} k \log_2(3/2)$, while preserving the full non-linear constraint for generic $v_2(k)$.

---

## 3. Dynamical & Automata-Theoretic Invariants

* **Carry Transducer on Continuous Ones:**
  The transition matrix on states $(a_{j-1}, c_j) \in \{00, 01, 10, 11\}$:
  $$A = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 1 & 0 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$
  satisfies $\det(A - \lambda I) = -\lambda (1 - \lambda)^2 (1 + \lambda)$. The spectral radius is $\rho(A) = 1$, yielding **zero topological entropy** $h_{\mathrm{top}}(\Sigma_A) = \log_2 \rho(A) = 0$. Combinatorially, the number of input preimages of length $m$ generating continuous ones grows at most polynomially ($N(m) \le C m$), rather than exponentially.
* **Transfer Operator on the $\beta = 3/2$ Shift:**
  The Perron–Frobenius operator $\mathcal{L}$ on $BV([0, 1])$ has essential spectral radius $r_{\mathrm{ess}}(\mathcal{L}) \le 2/3$. By Borel–Cantelli, the set of exceptional orbits has Lebesgue measure zero and Hausdorff dimension zero.

---

## 4. Analysis of Attempted Proof Extensions & Diophantine Barriers

1. **Frey Curve / Szpiro Bound:** Forces $\sigma_{\mathrm{fail}} \ge 5.26189$, whereas unconditional modular degree bounds yield only $\sigma \le 12$.
2. **Reciprocal Hermite–Padé Degeneracy:** The simultaneous Padé system for $(1-z)^{\pm 1/2}$ at $z = 1/9$ drops rank over $\mathbb{Q}$ ($3 f_1 - \frac{8}{3} f_2 = 0$), collapsing to Bennett's exponent $\lambda = 0.787$.
3. **Asymmetric Padé vs. Chebyshev Clearing Factor:** The clearing factor $(e^{1/3})^k \approx (1.3956)^k$ combined with height $(1.5)^k$ yields base $2.0934$, strictly exceeding the 2-adic divisibility base $2.0000$.
4. **Archimedean vs. 2-Adic Size Leakage:** The 2-adic modulus $2^{v_2(k)+3} \le 8k$ leaves $\asymp (1.5)^k / (8k) \to \infty$ candidate failure integers in $[1, q_k)$.
5. **Carry Transducer Suffix Inversion & Digit Complexity:** Inverting the scalar transducer $T_u$ across a run of ones $1^L$ ($L \approx 0.415 k$) yields $u B \equiv -(\gamma_{M-L} + 1) \pmod{2^L}$. We separate this analysis into two rigorous results:
   * **Proposition 6.1 (Qualitative Finiteness):** A subword run-length bound $L(u) \le c u$ with $c < \log_2(4/3) \approx 0.415037$ implies at most finitely many counterexamples.
   * **Theorem 6.2 (Conditional Effective Elimination):** An effective bound $L(u, M) < \lambda_{\mathrm{target}} M - \log_2 M$ for all windows $M \ge M_0$ with $M_0 \le 471{,}599{,}900$, combined with Kubina–Wunderlich (1990), establishes zero exceptions for all $k \ge 2$.

---

## 5. Artifact Directory Structure

* [`two_adic_proof.tex`](./two_adic_proof.tex) / [`two_adic_proof.pdf`](./two_adic_proof.pdf): Full 13-page publication manuscript with Table 1 benchmark comparison, complete proofs, Table of Contents, Appendix A with formal Lean 4 declarations, and references.
* [`formalization/TwoAdicDefectRigidity.lean`](./formalization/TwoAdicDefectRigidity.lean): 100% axiom-free Lean 4 formalization.
* [`two_adic_verification.py`](./two_adic_verification.py): Python test suite verifying all 7 theorems, carry partitions, simulated failure escapes, and automaton spectra up to $k = 100$.
