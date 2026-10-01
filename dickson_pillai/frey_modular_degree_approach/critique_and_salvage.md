# Ruthless Peer Review & Salvage Protocol: The Frey–Szpiro Framework for Dickson–Pillai

## 1. Fresh-Eyes Adversarial Critique of `frey_proof.tex`

We subject the Frey–Hellegouarch modular degree proof pipeline to an elite mathematical audit.

---

### Failure Mode 1: The Hoffstein–Lockhart Modular Degree Gap
- **The Claim in the Paper:**
  The manuscript asserts that the modular degree $\deg(\varphi_k)$ and symmetric square $L$-function $L(1, \operatorname{Sym}^2 f_k)$ produce an effective bound $\ln |\Delta(E_k)| \le 4.80 \ln N(E_k)$, which would contradict $\sigma_{\mathrm{fail}} \ge 5.26189$.
- **The Mathematical Reality:**
  Unconditionally, the sharpest known analytic bounds for symmetric square $L$-functions (Hoffstein–Lockhart 1994, Goldfeld–Hoffstein–Lieman 1994) give:
  $$c(\varepsilon) N_k^{-\varepsilon} \le L(1, \operatorname{Sym}^2 f_k) \le C(\varepsilon) (\ln N_k)^3$$
  When substituted into Faltings' height inequality:
  $$h_F(E_k) \le \frac{1}{2} \ln \deg(\varphi_k) + \frac{1}{2} \ln N_k \le \ln N_k + \frac{k}{2} \ln 3$$
  Multiplying by 12 yields:
  $$\ln |\Delta(E_k)| \le 12 \ln N(E_k) \implies \sigma(E_k) \le 12$$
- **The Defect:**
  The unconditional bound $\sigma(E_k) \le 12$ is **strictly weaker** than the failure threshold $\sigma_{\mathrm{fail}} \ge 5.26189$. An unconditional proof of $\sigma(E_k) < 5.26$ for all elliptic curves is currently equivalent to an explicit form of **Szpiro's Conjecture** ($s < 5.26$), which remains unproven.

---

### Failure Mode 2: Non-Uniformity of Period Integrals $\Omega_1(E_k)$
- The real period $\Omega(E_k) = \int_0^\infty \frac{dx}{\sqrt{x(x+a_k)(x+b_k)}} \approx \frac{\pi}{\sqrt{c_k}} = \frac{\pi}{3^{k/2}}$ decays exponentially.
- Because $\Omega(E_k) \to 0$, the modular degree $\deg(\varphi_k) \asymp \frac{\|f_k\|^2}{\Omega(E_k)^2} \asymp N_k \cdot 3^k$ is **exponentially large in $k$**.
- This exponential growth in $\deg(\varphi_k)$ introduces the exact $+ k \ln 3$ term that expands the Szpiro bound from $6$ to $12$.

---

## 2. The Salvage Architecture: Galois Level-Lowering (Mazur–Ribet on $E_k$)

To salvage the proof without needing the full, unproven Szpiro conjecture, we must exploit the **highly specialized prime-power structure** of the Dickson–Pillai partition:
$$r_k + 2^k q_k = 3^k$$
Notice that **two of the three terms ($2^k q_k$ and $3^k$) have massive prime-power valuations**:
$$v_2(B_k) = k, \quad v_3(C_k) = k$$

```
                           THE MAZUR-RIBET SALVAGE
                           
  E_k: y^2 = x(x - r_k)(x + 2^k q_k) ───> ρ_{E_k, p} : Gal(Qbar/Q) -> GL_2(F_p)
                                                   │
                                                   ▼  (v_2(Δ) = 2k, v_3(Δ) = 2k)
                                        Unramified at 2 and 3 for p | k
                                                   │
                                                   ▼  (Ribet's Level Lowering)
                                        g ∈ S_2(Γ_0(N_p)) with N_p = 2 rad(r_k)
```

### The Mechanism:
1. **Galois Representations mod $p$:**
   Let $p$ be a prime divisor of $k$ (or a prime $p \ge 5$). The mod-$p$ Galois representation is:
   $$\rho_{E_k, p} : \operatorname{Gal}(\overline{\mathbb{Q}}/\mathbb{Q}) \longrightarrow \operatorname{GL}_2(\mathbb{F}_p)$$
2. **Serre's Conductor & Ribet's Theorem:**
   Because $p \mid v_2(\Delta(E_k)) = 2k$ and $p \mid v_3(\Delta(E_k)) = 2k$, the representation $\rho_{E_k, p}$ is **unramified at 2 and 3** (except for the tame Serre inertia at $p$).
3. **The Stripped Conductor:**
   By Ken Ribet's Level-Lowering Theorem (1990), $\rho_{E_k, p}$ arises from a modular newform $g \in S_2(\Gamma_0(N_p))$ of reduced level:
   $$N_p = 2 \cdot \operatorname{rad}(r_k)$$
   **The exponential power $3^k$ and the quotient $q_k$ have been completely stripped from the modular level!**

4. **The Contradiction Threshold:**
   If $r_k + q_k > 2^k$, then $r_k$ is forced to satisfy the linear recurrence constraints imposed by the Hecke eigenvalues $a_\ell(g)$ of the fixed low-level space $S_2(\Gamma_0(2 \operatorname{rad}(r_k)))$.
   Because $\dim S_2(\Gamma_0(2 \operatorname{rad}(r_k))) \ll \operatorname{rad}(r_k)$, the modular forms cannot match the chaotic 2-adic expansion of $3^k \bmod 2^k$, yielding a finite bound on $k$.

---

## 3. Summary Assessment

| Approach | Status | Exact Obstacle |
| :--- | :--- | :--- |
| **Naive Szpiro / Modular Degree** | Needs $abc$ | Unconditional $L(1, \operatorname{Sym}^2)$ gives $\sigma \le 12$, missing $\sigma \le 5.26$ |
| **Ribet Level-Lowering on $(2^k, 3^k, r_k)$** | **Plausible & Active** | Strips $3^k, 2^k$ from level; bounds $r_k$ via Hecke eigenvalues |
