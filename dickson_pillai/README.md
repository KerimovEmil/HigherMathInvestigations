# The Dickson–Pillai Condition and Waring's Problem

## 1. Overview & Historical Background

In 1770, Edward Waring posed his celebrated problem: does there exist for every integer $k \ge 2$ a number $s$ such that every natural number is the sum of at most $s$ positive $k$-th powers? David Hilbert established the qualitative existence of $g(k)$ in 1909 (the Hilbert–Waring theorem).

In 1936, **Leonard Eugene Dickson** and **Subbayya Sivasankaranarayana Pillai** independently established that the exact value of $g(k)$ for all natural numbers is given by:

$$g(k) = 2^k + \lfloor (3/2)^k \rfloor - 2$$

provided that the following condition holds for all $k \ge 1$:

$$2^k \cdot \{(3/2)^k\} + \lfloor (3/2)^k \rfloor \le 2^k$$

This inequality is known as the **Dickson–Pillai Condition**.

---

## 2. Algebraic Formulations & Equivalences

Let $3^k = q_k 2^k + r_k$, with quotient $q_k = \lfloor (3/2)^k \rfloor$ and remainder $r_k = 3^k \bmod 2^k$ ($0 \le r_k < 2^k$). The fractional part is $\{(3/2)^k\} = r_k / 2^k$.

The condition can be rewritten in four strictly equivalent forms:

1. **Remainder & Quotient Form:**
   $$r_k + q_k \le 2^k$$

2. **Floor Form:**
   $$\lfloor (3/2)^k \rfloor \ge \frac{3^k - 2^k}{2^k - 1} = (3/2)^k \cdot \frac{1 - (2/3)^k}{1 - 2^{-k}}$$

3. **Ceiling Form:**
   $$\lceil (3/2)^k \rceil \ge \frac{3^k - 1}{2^k - 1}$$

4. **Fractional Part / Diophantine Distance Form:**
   $$1 - \{(3/2)^k\} \ge \frac{3^k - 2^k}{4^k - 2^k} = (3/4)^k \cdot \left(\frac{1 - (2/3)^k}{1 - 2^{-k}}\right) \approx (3/4)^k$$

---

## 3. Theoretical Landscape & Inherent Obstructions

* **Modular & 2-Adic Dynamics:** The residue $r_k = 3^k \bmod 2^k$ represents the low $k$ bits of $3^k$. Under the 2-adic multiplication map $x \mapsto 3x$ on $\mathbb{Z}_2$, the sequence behaves pseudo-randomly.
* **Mahler's Finiteness Theorem (1957):** Using Ridout's $p$-adic extension of Roth's theorem, Kurt Mahler proved that $1 - \{(3/2)^k\} > (3/4)^k$ holds for all but at most finitely many $k$. However, Roth's method is **ineffective** (provides no computable threshold $K_0$).
* **Baker's Effective Logarithmic Forms:** Effective Diophantine methods (Baker, Beukers 1981, Dubickas) yield lower bounds of the form $\|(3/2)^k\| > 2^{-c k}$. Currently, the best effective exponent is $c \approx 0.999$, whereas proving the conjecture requires $c \le \log_2(4/3) \approx 0.415037$.
* **Prime-Factor Leakage (Hermite–Padé Obstruction):** $r_k$ introduces prime factors outside any finite smooth base (e.g. $r_8 = 161 = 7 \times 23$), invalidating naive $S$-unit product formula closures.
* **Transcendence & Cusp Degeneration (Nesterenko Modular Obstruction):** For rational $q_k$, Eisenstein values $E_2(q_k), E_4(q_k), E_6(q_k)$ are algebraically independent and transcendental (Nesterenko 1996), so $3^{kL} \theta^m \Phi(q_k) \notin \mathbb{Z}$. Rigorous Philippon elimination resultants and Siegel's lemma height penalties restore the classical Baker barrier.
* **Computational Boundary:** In 1990, Kubina & Wunderlich verified that **zero exceptions** exist for all $k \le 471{,}600{,}000$.

---

## 4. Connections to Related Open & Closed Problems

| Problem | Key Formulation | Nature | Status |
| :--- | :--- | :--- | :--- |
| **Dickson–Pillai Condition** | $1 - \{(3/2)^k\} \ge (3/4)^k$ | Powers of $3/2$ near 1 | Finitely many exceptions (Mahler 1957); **0 exceptions** conjectured |
| **Mahler's $3/2$ Problem** | $\{\xi (3/2)^n\} < 1/2$ | Bounded fractional parts ($Z$-numbers) | **Open** |
| **Equidistribution of $(3/2)^n$** | $\{(3/2)^n\} \bmod 1$ | Dense / uniform distribution | **Open** |
| **Collatz ($3x+1$) Problem** | $3^k \bmod 2^k$ | 2-adic bit dynamics | **Open** |
| **Pisot–Vijayaraghavan Numbers** | $\mathrm{dist}(\alpha^n, \mathbb{Z}) \to 0$ | Exponential convergence to integers | **Solved** (Pisot 1938) |
| **$abc$ Conjecture** | $C \le K(\varepsilon) \mathrm{rad}(ABC)^{1+\varepsilon}$ | Prime power additive constraints | **Open** |

---

## 5. Directory Structure & Research Organization

The investigation is organized into four clean modular subdirectories:

### 📁 `manuscript_and_formalization/` (Unified Theoretical Paper & Formalization)
* [`paper.tex`](./manuscript_and_formalization/paper.tex) / [`paper.pdf`](./manuscript_and_formalization/paper.pdf): Comprehensive peer-reviewed research paper detailing both Hermite–Padé and Automorphic frameworks, along with their Diophantine obstructions.
* [`formalization/DicksonPillai.lean`](./manuscript_and_formalization/formalization/DicksonPillai.lean): Axiom-free Lean 4 formalization proving base cases, prime leakage counterexamples, $\mathrm{SL}_2(\mathbb{Z})$ modular transformations, Ramanujan derivations, and Diophantine exponent hierarchies.
* [`5adic_elimination_and_k0.md`](./manuscript_and_formalization/5adic_elimination_and_k0.md): Detailed mathematical analysis of the 5-adic elimination bounds and prime leakage.
* [`explicit_cutoff_k0.py`](./manuscript_and_formalization/explicit_cutoff_k0.py): Non-asymptotic constant scripts.
* [`hermite_pade_systematic.py`](./manuscript_and_formalization/hermite_pade_systematic.py): Multi-dimensional simultaneous Hermite–Padé engine.
* [`pade_hypergeometric.py`](./manuscript_and_formalization/pade_hypergeometric.py): Padé approximant engine for hypergeometric forms.

### 📁 `nesterenko_modular_approach/` (Automorphic & Eisenstein Differential Algebra)
* [`modular_proof.tex`](./nesterenko_modular_approach/modular_proof.tex) / [`modular_proof.pdf`](./nesterenko_modular_approach/modular_proof.pdf): Dedicated manuscript examining Nesterenko's differential algebra $(\mathbb{Q}[q, E_2, E_4, E_6], \theta)$ and analyzing the three Diophantine failure modes.
* [`modular_forms.py`](./nesterenko_modular_approach/modular_forms.py): High-precision evaluation of Eisenstein series $E_2, E_4, E_6$ and Ramanujan derivations to $> 100$ digits.
* [`auxiliary_modular_polynomial.py`](./nesterenko_modular_approach/auxiliary_modular_polynomial.py): Monomial basis and kernel solver for auxiliary modular forms.
* [`cusp_inversion_bound.py`](./nesterenko_modular_approach/cusp_inversion_bound.py): Quasi-modular $S$-inversion asymptotics near $\tau \to 0$.
* [`stress_test.py`](./nesterenko_modular_approach/stress_test.py): Rigorous validation suite for modular differential algebra.
* [`uniform_zero_lemma.py`](./nesterenko_modular_approach/uniform_zero_lemma.py): Multiplicity and zero-order calculations.

### 📁 `frey_modular_degree_approach/` (Frey–Hellegouarch Curves & Szpiro Framework)
* [`frey_proof.tex`](./frey_modular_degree_approach/frey_proof.tex) / [`frey_proof.pdf`](./frey_modular_degree_approach/frey_proof.pdf): Research manuscript detailing the Frey curve $E_k : y^2 = x(x - r_k)(x + 2^k q_k)$ and Szpiro asymptotic thresholds.
* [`formalization/FreyDicksonPillai.lean`](./frey_modular_degree_approach/formalization/FreyDicksonPillai.lean): Verified Lean 4 formalization of Frey curve Weierstrass models, full 2-torsion, and Szpiro threshold inequalities.
* [`critique_and_salvage.md`](./frey_modular_degree_approach/critique_and_salvage.md): Adversarial critique of the modular degree gap and the Mazur–Ribet level-lowering salvage protocol.

### 📁 `two_adic_rigidity_approach/` (2-Adic Valuation Rigidity & Carry Drift)
* [`two_adic_proof.tex`](./two_adic_rigidity_approach/two_adic_proof.tex) / [`two_adic_proof.pdf`](./two_adic_rigidity_approach/two_adic_proof.pdf): Comprehensive research paper detailing exact LTE valuation formulas, carry drift dynamics, and the Archimedean vs.\ 2-adic barrier.
* [`formalization/TwoAdicDefectRigidity.lean`](./two_adic_rigidity_approach/formalization/TwoAdicDefectRigidity.lean): 100% axiom-free, sorry-free Lean 4 formalization proving division identities, multiplier spectrum bounds, parity valuations, and candidate leakage.
* [`two_adic_verification.py`](./two_adic_rigidity_approach/two_adic_verification.py): High-precision Python verification script auditing valuations and multiplier spectra.
* [`README.md`](./two_adic_rigidity_approach/README.md): Modular framework documentation.

### 📁 `numerical_validation/` (Streaming Engines & High-Throughput Verification)
* [`parallel_verifier.cpp`](./numerical_validation/parallel_verifier.cpp): 20-thread parallel C++ engine with unrolled 128-bit limb arithmetic.
* [`verifier.cpp`](./numerical_validation/verifier.cpp): Single-threaded reference verifier with rolling state hashing.
* [`dickson_pillai.py`](./numerical_validation/dickson_pillai.py): Arbitrary-precision integer verifier and statistical analyzer.
* [`verification_records/`](./numerical_validation/verification_records/):
  * [`checkpoints_25M.json`](./numerical_validation/verification_records/checkpoints_25M.json) / [`summary_25M.md`](./numerical_validation/verification_records/summary_25M.md): Certified parallel ledger for $k = 1$ to $25{,}000{,}000$.
  * [`checkpoints_10M.json`](./numerical_validation/verification_records/checkpoints_10M.json) / [`summary_10M.md`](./numerical_validation/verification_records/summary_10M.md): Checkpoints for $k \le 10{,}000{,}000$.
  * [`checkpoints_5M.json`](./numerical_validation/verification_records/checkpoints_5M.json) / [`summary_5M.md`](./numerical_validation/verification_records/summary_5M.md): Checkpoints for $k \le 5{,}000{,}000$.

### 📁 `miscellaneous/` (Plots & Supporting Visualizations)
* [`plots/dickson_pillai_analysis.png`](./miscellaneous/plots/dickson_pillai_analysis.png): Safety margin growth and danger envelope proximity.
* [`plots/pade_exponent_analysis.png`](./miscellaneous/plots/pade_exponent_analysis.png): Padé error decay and effective exponent trajectories.
* [`plots/hermite_pade_candidate_analysis.png`](./miscellaneous/plots/hermite_pade_candidate_analysis.png): Comparative decay rates across algebraic candidate points.
* [`plots/explicit_k0_crossover.png`](./miscellaneous/plots/explicit_k0_crossover.png): Logarithmic lower bound vs barrier crossover.
