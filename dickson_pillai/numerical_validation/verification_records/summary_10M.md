# Parallel Verification Run Summary: $k = 1$ to $10{,}000{,}000$

- **Date of Execution:** 2026-09-30
- **Engine:** C++ Parallel Multi-Core 2-Adic Streaming Verifier (`parallel_verifier.cpp`)
- **Hardware:** AMD Ryzen AI 9 365 (10 physical Zen 5 cores, 20 logical threads)
- **Range Checked:** $k \in [1, 10{,}000{,}000]$
- **Total Exceptions / Violations:** **0**
- **Verification Status:** **100% Valid (Passed)**
- **Total Execution Time:** 170.376 seconds (2 min 50 sec)
- **Aggregate Throughput:** **58,693 steps/second**
- **Chunks Processed:** 20 chunks (500,000 steps per chunk)

---

## Order Statistics (Top Extreme Near-Misses)

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
