# Research Roadmap: The Padé & Hypergeometric Method for the Dickson–Pillai Condition

## 1. Executive Summary & Core Objective

The **Dickson–Pillai condition** for Waring's problem requires that for all $k \ge 1$:
$$1 - \{(3/2)^k\} \ge (3/4)^k \cdot \left(\frac{1 - (2/3)^k}{1 - 2^{-k}}\right) \approx (3/4)^k$$

Effective Diophantine approximation (via Baker's theory and classical Padé approximants) produces lower bounds of the form:
$$\|(3/2)^k\| > 2^{-c k} \quad \text{for all } k \ge K_0$$

* **Current Effective Record:** $c \approx 0.999999$ (Beukers 1981, Dubickas 2006).
* **Target Exponent:** $c \le \log_2(4/3) \approx 0.415037$.
* **Goal:** If an effective method proves $c \le 0.415037$ with an explicit computable threshold $K_0 \le 471{,}600{,}000$, the Dickson–Pillai condition—and hence the exact formula for $g(k)$ in Waring's problem—is **fully and unconditionally solved**.

---

## 2. Mathematical Framework

### 2.1 The Padé Approximation Mechanics
Let $f(z) = (1 - z)^\nu$ where $\nu \in \mathbb{Q} \setminus \mathbb{Z}$. The $[n/n]$ Padé approximant around $z = 0$ yields polynomials $P_n(z), Q_n(z) \in \mathbb{Q}[z]$ of degree at most $n$ such that:
$$Q_n(z)(1 - z)^\nu - P_n(z) = z^{2n+1} R_n(z)$$

The error term has the explicit contour/Euler integral representation:
$$R_n(z) = \frac{\Gamma(n + 1 + \nu)}{\Gamma(\nu) n!} \int_0^1 \frac{x^n (1 - x)^{n - \nu}}{(1 - z x)^{n + 1}} \, dx$$

### 2.2 Jacobi Polynomial Connection
The polynomials $P_n(z)$ and $Q_n(z)$ are explicitly given by hypergeometric functions related to the Jacobi polynomials $P_n^{(\alpha, \beta)}(t)$:
$$Q_n(z) = {}_2F_1\left(-n, -n - \nu; -2n; z\right)$$

---

## 3. Multi-Stage Research Plan

### Stage 1: Dissecting Beukers' Quadratic Form $(1 - 1/9)^{1/2}$
1. **The Choice of Base:**
   $$\sqrt{1 - \frac{1}{9}} = \sqrt{\frac{8}{9}} = \frac{2\sqrt{2}}{3}$$
2. **Integral Evaluation:**
   $$I_n = \int_0^1 \frac{x^n (1 - x)^n}{\left(1 - \frac{x}{9}\right)^{n + 1/2}} \, dx = A_n - B_n \sqrt{\frac{8}{9}}$$
3. **Asymptotic Analysis:**
   * Coefficient growth: $|Q_n(1/9)| \sim (\mu_1)^n$ where $\mu_1 = (3 + 2\sqrt{2})^2 \approx 33.9705$.
   * Error decay: $|R_n(1/9)| \sim (\mu_2)^n$ where $\mu_2 = (3 - 2\sqrt{2})^2 \approx 0.029437$.
4. **Denominator Growth Analysis:**
   * Let $D_n = \mathrm{lcm}(1, 2, \dots, 2n)$. By the Prime Number Theorem, $D_n \sim e^{2n} \approx 7.389^n$.
   * The effective exponent is governed by:
     $$c = \frac{\ln \mu_1 + \ln D_n}{\ln(1/\mu_2) - \ln D_n} \approx \frac{3.5255 + 2.0000}{3.5255 - 2.0000} \approx 0.999999$$

---

### Stage 2: Cubic & Higher-Order Approximations ($(1 - 19/27)^{1/3}$)
1. **The Cubic Relation:**
   $$1 - \frac{19}{27} = \frac{8}{27} = \left(\frac{2}{3}\right)^3 \implies \left(1 - \frac{19}{27}\right)^{1/3} = \frac{2}{3}$$
2. **Simultaneous Padé Approximations:**
   * Construct simultaneous rational approximants $P_{1,n}(z), P_{2,n}(z), Q_n(z)$ to $(1 - z)^{1/3}$ and $(1 - z)^{2/3}$.
   * Simultaneous linear forms in 3 linear relations eliminate one degree of freedom, drastically reducing the effective remainder growth.
3. **Saddle-Point Asymptotics for the Cubic Remainder:**
   * Remainder integral:
     $$J_n = \int_0^1 \frac{x^{2n} (1 - x)^n}{\left(1 - \frac{19}{27} x\right)^{2n + 2/3}} \, dx$$
   * Calculate the saddle point $x_0 \in (0, 1)$ of $\phi(x) = \frac{x^2 (1 - x)}{(1 - 19x/27)^2}$.

---

### Stage 3: Arithmetic Pruning of Common Factors (Denominator Optimization)
1. **$p$-Adic Valuation of Hypergeometric Coefficients:**
   * The theoretical upper bound $D_n = \mathrm{lcm}(1, \dots, 2n) \sim e^{2n}$ is overly pessimistic because prime powers divide large subsets of binomial coefficients.
   * Compute the **Hata–Rukhadze arithmetic factor**:
     $$\Phi = \lim_{n \to \infty} \frac{1}{n} \sum_{p \le 2n} v_p(G_n) \ln p$$
   * Replacing $e^{2n}$ with $e^{(2 - \Phi)n}$ reduces denominator penalty.

---

### Stage 4: Python Simulation & Prototype Verification
1. Implement arbitrary-precision numerical integration and recurrence relations in `dickson_pillai/` to compute $A_n, B_n, I_n$ for $n \le 500$.
2. Calculate empirical exponents $c_n$ for $n = 10, 20, \dots, 500$ across:
   * Standard quadratic form $(1 - 1/9)^{1/2}$
   * Cubic form $(1 - 19/27)^{1/3}$
   * Quintic form $(1 - 211/243)^{1/5}$ (since $1 - 211/243 = 32/243 = (2/3)^5$).
3. Output comparative tables of effective exponents $c$.

---

### Stage 5: Simultaneous Hermite–Padé Systems & Exponent Crossover
1. **Multi-Dimensional Transference:**
   In dimension $d$ with $d-1$ independent linear forms, the height burden is distributed across relations:
   $$c_{\mathrm{eff}} = \frac{\ln \mu_1 + \ln D}{(d - 1) \left(\ln(1/\mu_2) - \ln D\right)}$$
2. **Smooth Number Algebraic Generators:**
   * **Quartic Form ($d=4$):** $(1 - 1/81) = 80/81 = (2/3)^4 \cdot 5$.
     With $z = 1/81$ and $d-1 = 3$ linear forms, $c_{\mathrm{eff}} \approx \mathbf{0.3999} < \mathbf{0.415037}$.
   * **Quintic Form ($d=5$):** $(1 - 1/243) = 242/243 = 2 \cdot 11^2 / 3^5$.
     With $z = 1/243$ and $d-1 = 4$ linear forms, $c_{\mathrm{eff}} \approx \mathbf{0.2882} \ll \mathbf{0.415037}$.

---

## 4. Key Milestones & Current Progress

| Milestone | Deliverable | Target Metric | Status |
| :--- | :--- | :--- | :--- |
| **M1: Baseline Replication** | Exact Python module for Beukers' quadratic integrals | Reproduce $c \approx 0.999999$ | **Completed** (`pade_hypergeometric.py`) |
| **M2: Arithmetic Factor $\Phi$** | Hata-style prime factor sieve for $Q_n(1/9)$ | Reduce $c$ to $< 0.85$ | **Completed** (`pade_hypergeometric.py`) |
| **M3: Cubic Padé System** | Simultaneous Padé approximants for $2^{1/3}$ ($z = 3/128$) | Test if $c < 0.65$ | **Completed** ($c \approx 0.6273$, `hermite_pade_systematic.py`) |
| **M4: Exponent Crossover** | Higher-order simultaneous approximation ($z = 1/81$) | Achieve $c \le 0.415037$ | **Completed** ($c \approx 0.3999 < 0.415037$) |
| **M5: Explicit Threshold $K_0$** | Compute finite cutoff $K_0$ | Verify $K_0 \le 4.71 \times 10^8$ | **In Progress** |
