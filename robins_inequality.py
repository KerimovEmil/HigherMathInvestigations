"""
https://en.wikipedia.org/wiki/Divisor_function
"""
from math import exp, log
from divisor import divisors
import matplotlib.pyplot as plt


def harmonic(n: int) -> float:
    """Computes the nth harmonic number"""
    h = 0.0
    for i in range(1, n+1):
        h += 1.0/i
    return h


def robin_ineq(n):
    lhs = divisors(1, n)
    h_n = harmonic(n)
    rhs = h_n + exp(h_n) * log(h_n)
    return rhs - lhs


if __name__ == '__main__':
    n_range = range(5041, 10000)
    ls = [robin_ineq(n) for n in n_range]

    plt.scatter(n_range, ls, s=5, label='h_n + e^{h_n}*ln(h_n) - sigma(n)')
    plt.title("Robin's inequality")
    plt.xlabel("n")
    plt.ylabel("h_n + e^{h_n}*ln(h_n) - sigma(n)")
    plt.show()
