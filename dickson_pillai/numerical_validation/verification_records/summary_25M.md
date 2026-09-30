# Parallel Verification Run Summary: $k = 1$ to $25{,}000{,}000$

- **Date of Execution:** 2026-09-30
- **Engine:** C++ Parallel Multi-Core 2-Adic Streaming Verifier (`parallel_verifier.cpp`)
- **Hardware:** AMD Ryzen AI 9 365 (10 physical Zen 5 cores, 20 logical threads)
- **Range Checked:** $k \in [1, 25{,}000{,}000]$
- **Total Exceptions / Violations:** **0**
- **Verification Status:** **100% Valid (Passed)**
- **Total Execution Time:** 1710.354 seconds (28 min 30 sec)
- **Aggregate Throughput:** **14,616 steps/second** (at $k \approx 2.5 \times 10^7$)
- **Chunks Processed:** 50 chunks (500,000 steps per chunk)
- **Final Memory Size of $3^k$:** 619,130 64-bit words (~39.6 million bits)

---

## 1. Mathematical Dynamics: Why Global Near-Misses Cluster at Small $k$

The **Safety Ratio** is defined as:
$$\text{Safety Ratio}(k) = \frac{1 - \{(3/2)^k\}}{(3/4)^k} \approx (1 - \theta_k) \cdot \left(\frac{4}{3}\right)^k$$

Because $(4/3)^k \approx 2^{0.415 k}$ grows exponentially with $k$:
* At $k = 14$: $(4/3)^{14} \approx 56.12 \implies \text{Safety Ratio} = 3.9701$.
* At $k = 100$: $(4/3)^{100} \approx 3.1 \times 10^{12} \implies$ requires $1 - \theta_{100} < 10^{-12}$ to beat $k = 14$.
* At $k = 25{,}000{,}000$: $(4/3)^{25,000,000} \approx 10^{3,123,000} \implies$ requires $\theta_k$ within $10^{-3,123,000}$ of $1$.

Because fractional parts $\{(3/2)^k\}$ are uniformly distributed in $[0, 1)$, the probability of beating the small-$k$ safety ratios vanishes exponentially. Thus, the global order statistics list **theoretically must be dominated by small $k \le 16$**, confirming the soundness of the distribution.

---

## 2. Order Statistics (Global Extreme Near-Misses)

| Rank | $k$ | Remainder Fraction $\{(3/2)^k\}$ | Safety Ratio $\frac{1 - \theta_k}{(3/4)^k}$ | $\log_2(\text{Safety Ratio})$ |
| :--- | :--- | :--- | :--- | :--- |
| **1** | 2 | 0.250000 | 1.3333 | 0.4150 |
| **2** | 3 | 0.375000 | 1.4815 | 0.5670 |
| **3** | 5 | 0.593750 | 1.7119 | 0.7756 |
| **4** | 4 | 0.062500 | 2.9630 | 1.5670 |
| **5** | 6 | 0.390625 | 3.4239 | 1.7756 |
| **6** | 8 | 0.628906 | 3.7068 | 1.8902 |
| **7** | 14 | 0.929260 | 3.9701 | 1.9892 |
| **8** | 10 | 0.665039 | 5.9481 | 2.5724 |
| **9** | 7 | 0.085938 | 6.8477 | 2.7756 |
| **10** | 9 | 0.443359 | 7.4135 | 2.8902 |
| **11** | 15 | 0.893890 | 7.9403 | 2.9892 |
| **12** | 12 | 0.746338 | 8.0079 | 3.0014 |
| **13** | 11 | 0.497559 | 11.8963 | 3.5724 |
| **14** | 16 | 0.840836 | 15.8806 | 3.9892 |
| **15** | 13 | 0.619507 | 16.0159 | 4.0014 |

---

## 3. Scale Progression Comparison

| Verification Dataset | Range Verified | Execution Time | Limb Count (64-bit) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **5 Million** | $k \le 5 \times 10^6$ | 278.36 s | 123,826 | Verified |
| **10 Million** | $k \le 10 \times 10^6$ | 170.38 s | 247,652 | Verified |
| **25 Million** | $k \le 25 \times 10^6$ | 1710.35 s | 619,130 | Verified |
