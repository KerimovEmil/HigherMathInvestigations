import random
import math
from tqdm import tqdm

"""
PROOF: Expected Value of the Stopping Increment Y_N

Problem:
Sample Y_i ~ U(0, 1) independently. 
Stop at N = min{ n : sum(Y_1...Y_n) >= 1 }.
Find E[Y_N].

-------------------------------------------------------------------------------
1. EXPECTED STOPPING TIME (E[N])
The probability that the sum of n uniform variables is less than 1 is the 
volume of an n-dimensional simplex: P(S_n < 1) = 1/n!.
Since N > n is equivalent to S_n < 1:
    E[N] = sum_{n=0 to inf} P(N > n) = sum_{n=0 to inf} 1/n! = e.

2. WALD'S IDENTITY FOR THE TOTAL SUM (E[S_N])
Wald's Identity states E[S_N] = E[N] * E[Y_i].
    E[S_N] = e * (1/2) = e/2.

3. EXPECTED SUM BEFORE STOPPING (E[S_{N-1}])
For n >= 2, the density of S_{n-1} at s (for s < 1) is f(s) = s^(n-2)/(n-2)!.
The probability of stopping at step n given S_{n-1} = s is P(s + Y_n >= 1) = s.
    E[S_{N-1}] = sum_{n=2 to inf} integral_0^1 [ s * (s) * s^(n-2)/(n-2)! ] ds
               = sum_{n=2 to inf} integral_0^1 [ s^n / (n-2)! ] ds
               = sum_{n=2 to inf} [ 1 / ((n+1)(n-2)!) ]

Let k = n - 2:
    E[S_{N-1}] = sum_{k=0 to inf} [ 1 / ((k+3)k!) ]
Using the power series integral of x^2 * e^x from 0 to 1:
    E[S_{N-1}] = e - 2.

4. THE LAST INCREMENT (E[Y_N])
Since S_N = S_{N-1} + Y_N:
    E[Y_N] = E[S_N] - E[S_{N-1}]
    E[Y_N] = (e/2) - (e - 2)
    E[Y_N] = 2 - e/2 ≈ 0.640859
-------------------------------------------------------------------------------
"""

NUM_SIMULATIONS = 100_000
theoretical = 2 - math.e / 2

last_values = []
running_avg = 0.0

pbar = tqdm(range(NUM_SIMULATIONS), desc="Simulating")

for i in pbar:
    total = 0.0
    while True:
        y = random.uniform(0, 1)
        total += y
        if total >= 1:
            last_values.append(y)
            break

    running_avg = sum(last_values) / (i + 1)
    pbar.set_postfix(avg=f"{running_avg:.7f}", theory=f"{theoretical:.7f}", diff=f"{abs(running_avg - theoretical):.7f}")

pbar.close()

print()
print(f"Simulated average of y_N : {running_avg:.7f}")
print(f"Theoretical (2 - e/2)    : {theoretical:.7f}")
print(f"Difference               : {abs(running_avg - theoretical):.7f}")

#