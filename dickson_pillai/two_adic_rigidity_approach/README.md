# 2-Adic Defect Rigidity & Carry Drift Framework

## Overview

This directory contains the theoretical manuscript, computational verification suite, and verified Lean 4 formalization of the **$2$-adic defect rigidity and discrete carry drift** framework for the Dickson–Pillai condition in Waring's problem:

$$r_k + q_k \le 2^k \iff D_k \ge q_k, \quad \text{where } 3^k = q_k 2^k + r_k \text{ and } D_k := 2^k - r_k$$

---

## Key Proven Results

1. **Even Index 2-Adic Rigidity ($k = 2m$):**
   Using $3^k - 1 = 9^m - 1 = 8 V_m$ and the Lifting the Exponent Lemma (LTE), for all even $k \ge 6$:
   $$v_2(D_k + 1) = v_2(k) + 2$$

2. **Odd Index 2-Adic Rigidity ($k \equiv 3 \pmod 4$):**
   For all $k \equiv 3 \pmod 4$ ($k \ge 5$):
   $$D_k \equiv 5 \pmod{16} \iff v_2(D_k + 3) = 3$$

3. **Odd Index Branch ($k \equiv 1 \pmod 4$):**
   For all $k \equiv 1 \pmod 4$ ($k - v_2(k-1) \ge 3$):
   $$v_2(D_k + 3) = v_2(k - 1) + 2$$

4. **Residue Exclusions & Small Defect Impossibility:**
   For all $k \ge 3$:
   $$D_k \bmod 8 \in \{5, 7\}$$
   Consequently, the small defect equations $D_k = 1$ and $D_k = 3$ have **zero solutions** for all $k \ge 3$. Furthermore, $D_k$ is strictly odd for all $k \ge 1$.

5. **Discrete Carry Drift Recurrence & Multiplier Spectrum:**
   For all $k \ge 1$:
   $$3 D_k - D_{k+1} = C_k 2^k, \quad C_k := 3 q_k - 2 q_{k+1} + 1 \in \{-1, 0, 1, 2\}$$
   The naive conjecture $C_k \in \{-1, 0, 1\}$ is formally disproved by the counterexamples $C_4 = C_7 = 2$.

6. **The Archimedean vs. 2-Adic Barrier:**
   The $2$-adic modulus is bounded by $2^{v_2(k)+3} \le 8k$ (linear/logarithmic scale), while the failure threshold $q_k \approx (1.5)^k$ grows exponentially. Within the danger window $[1, q_k)$, there exist:
   $$\asymp \frac{(1.5)^k}{8k} \to \infty$$
   candidate integers satisfying the exact $2$-adic rigidity equations. Consequently, purely non-Archimedean $2$-adic methods cannot establish the Archimedean inequality $D_k \ge q_k$.

---

## Directory Structure

* [`two_adic_proof.tex`](./two_adic_proof.tex) / [`two_adic_proof.pdf`](./two_adic_proof.pdf): Full research paper documenting all proofs, LTE derivations, carry recurrence analysis, and the Archimedean barrier.
* [`formalization/TwoAdicDefectRigidity.lean`](./formalization/TwoAdicDefectRigidity.lean): 100% axiom-free, sorry-free Lean 4 formalization.
* [`two_adic_verification.py`](./two_adic_verification.py): High-precision Python verification script validating all identities and multiplier values up to $k = 100$.
