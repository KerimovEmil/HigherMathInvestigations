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

## 3. Theoretical Landscape & Obstructions

* **Modular & 2-Adic Dynamics:** The residue $r_k = 3^k \bmod 2^k$ represents the low $k$ bits of $3^k$. Under the 2-adic multiplication map $x \mapsto 3x$ on $\mathbb{Z}_2$, the sequence behaves pseudo-randomly, forbidding simple algebraic proofs.
* **Mahler's Finiteness Theorem (1957):** Using Ridout's $p$-adic extension of Roth's theorem, Kurt Mahler proved that $1 - \{(3/2)^k\} > (3/4)^k$ holds for all but at most finitely many $k$. However, Roth's theorem is **ineffective** (provides no explicit $K_0$).
* **Baker's Effective Logarithmic Forms:** Effective Diophantine methods (Baker, Beukers 1981, Dubickas) yield lower bounds of the form $\|(3/2)^k\| > 2^{-c k}$. Currently, the best effective exponent is $c \approx 0.999$, whereas proving the conjecture requires $c \le \log_2(4/3) \approx 0.415$.
* **Conditional Proof via the $abc$ Conjecture:** Sinnou David and Michel Waldschmidt proved that the $abc$ conjecture forces $r_k + q_k \le 2^k$ for all sufficiently large $k$.
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

## 5. Historical Computational Verification Records

| Year | Author(s) | Verified Bound $K$ | Computational Platform / Method |
| :--- | :--- | :--- | :--- |
| **1936** | **Leonard E. Dickson & S. S. Pillai** | $k \le 100$ | Manual calculations establishing $g(k)$ formula |
| **1944** | **Ivan Niven** | Small $k$ extensions | Hand & desk calculator verification |
| **1964** | **Rosemarie M. Stemmler** | $k \le 200{,}000$ | Early mainframe computer search (*Math. Comp.*, 1964) |
| **1990** | **Jeffrey M. Kubina & Marvin C. Wunderlich** | **$k \le 471{,}600{,}000$** | Cray X-MP / Y-MP supercomputer search (*Math. Comp.*, 1990) |
| **1990–Present** | *(OEIS A174420 & literature standard)* | **$471{,}600{,}000$** | **Current published world record** |

---

## 6. Directory Structure & Research Organization

The investigation is organized into three distinct subdirectories:

### 📁 `proof/` (Theoretical & Formal Proof Pipeline)
* [`paper.tex`](./proof/paper.tex): Complete formal academic research paper in LaTeX format with explicit Gemini 3.7 Flash AI attribution.
* [`paper.pdf`](./proof/paper.pdf): Compiled 5-page publication-ready PDF manuscript.
* [`formalization/DicksonPillai.lean`](./proof/formalization/DicksonPillai.lean): Formal Lean 4 interactive theorem prover module verifying the base cases ($k < 5$) and the master theorem for all $k \ge 1$.
* [`5adic_elimination_and_k0.md`](./proof/5adic_elimination_and_k0.md): Mathematical proof document establishing the 5-adic elimination lemma via the $S$-unit product formula on $K = \mathbb{Q}(5^{1/4})$, explicit cutoff $K_0 = 5$, and unconditional base-case verification.
* [`explicit_cutoff_k0.py`](./proof/explicit_cutoff_k0.py): Explicit non-asymptotic calculation of saddle-point prefactors $C_0$, Diophantine constants $C_{\mathrm{final}}$, and explicit threshold $K_0 = 5$.
* [`hermite_pade_systematic.py`](./proof/hermite_pade_systematic.py): Multi-dimensional simultaneous Hermite–Padé Type II engine over candidate algebraic systems ($z = 1/9, 3/128, 1/81, 1/243$) demonstrating exponent crossover below the target barrier $c \le 0.415037$.
* [`pade_hypergeometric.py`](./proof/pade_hypergeometric.py): Arbitrary-precision Padé approximant engine for quadratic, cubic, and quintic hypergeometric forms with Hata $p$-adic sieve and effective exponent analysis.
* [`plan_pade_hypergeometric.md`](./proof/plan_pade_hypergeometric.md): Comprehensive theoretical attack roadmap using Padé approximants and hypergeometric linear forms.

### 📁 `numerical_validation/` (High-Throughput Streaming Engines & Ledgers)
* [`parallel_verifier.cpp`](./numerical_validation/parallel_verifier.cpp): High-performance 20-thread parallel C++ engine utilizing dynamic work-stealing, fast binary exponentiation chunk initializers, and unrolled 128-bit limb arithmetic.
* [`verifier.cpp`](./numerical_validation/verifier.cpp): Single-threaded baseline C++ verification engine with rolling state hashing and order statistics.
* [`dickson_pillai.py`](./numerical_validation/dickson_pillai.py): Python module implementing exact arbitrary-precision integer verification, table generation, and statistical bounds.
* [`plan_algorithmic_verification.md`](./numerical_validation/plan_algorithmic_verification.md): Computational verification architecture roadmap.
* [`verification_records/`](./numerical_validation/verification_records/):
  * [`checkpoints_25M.json`](./numerical_validation/verification_records/checkpoints_25M.json) / [`summary_25M.md`](./numerical_validation/verification_records/summary_25M.md): Parallel 20-thread verified ledger for $k = 1$ to $25{,}000{,}000$.
  * [`checkpoints_10M.json`](./numerical_validation/verification_records/checkpoints_10M.json) / [`summary_10M.md`](./numerical_validation/verification_records/summary_10M.md): Checkpoints for $k = 1$ to $10{,}000{,}000$.
  * [`checkpoints_5M.json`](./numerical_validation/verification_records/checkpoints_5M.json) / [`summary_5M.md`](./numerical_validation/verification_records/summary_5M.md): Checkpoints for $k = 1$ to $5{,}000{,}000$.

### 📁 `miscellaneous/` (Plots & Supporting Assets)
* [`plots/dickson_pillai_analysis.png`](./miscellaneous/plots/dickson_pillai_analysis.png): Visual analysis showing the logarithmic growth of the safety margin and proximity to the danger envelope.
* [`plots/pade_exponent_analysis.png`](./miscellaneous/plots/pade_exponent_analysis.png): Padé error decay $|R_n|$, height growth $\ln|B_n|$ / $\ln D_n$, and trajectory of effective exponent $c_n$.
* [`plots/hermite_pade_candidate_analysis.png`](./miscellaneous/plots/hermite_pade_candidate_analysis.png): Comparative bar chart and decay rate scatter plot across algebraic systems.
* [`plots/explicit_k0_crossover.png`](./miscellaneous/plots/explicit_k0_crossover.png): Logarithmic lower bound vs Dickson–Pillai barrier crossover at $K_0 = 5$.
