"""
phi(x) = sum_{p^n<=x} log(p) = sum_{k<=x} \Lambda(k)

phi(x) = x - 1/2 ln(1-x^-2) - ln(2pi) + 4*sqrt(x)* sum_{r} (cos(rln(x)) + 2r*sin(rln(x))) /(1+4r^2)
where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true)
"""

import numpy as np
import matplotlib.pyplot as plt
import primesieve
from math import log, sin, cos, pi


def von_mangoldt(n: int, primes_set: set[int]) -> float:
    """Compute the von Mangoldt function for an integer n using a set of primes."""
    if n in primes_set:
        return log(n)

    for prime in primes_set:
        power = prime * prime
        while power <= n:
            if power == n:
                return log(prime)
            power *= prime
    return 0


def chebyshev_psi(x):
    """Compute the Chebyshev psi function up to x."""
    primes = primesieve.primes(x)
    primes_set = set(primes)
    psi_values = np.zeros(x + 1)

    for n in range(1, x + 1):
        psi_values[n] = psi_values[n - 1] + von_mangoldt(n, primes_set)

    return psi_values


def plot_chebyshev_psi(x):
    """Plot the Chebyshev psi function up to x."""
    psi_values = chebyshev_psi(x)
    plt.figure(figsize=(10, 6))
    plt.plot(range(x + 1), psi_values, label=r'$\psi(x)$')
    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.title('Chebyshev $\psi$ Prime Counting Function')
    plt.legend()
    plt.grid(True)
    plt.show()


def chebyshev_psi_rz_approximation(x, ls_rz_zero):
    """
    Compute the Chebyshev psi function up to x, using
    phi(x) = x - 1/2 ln(1-x^-2) - ln(2pi) + 4*sqrt(x)* sum_{r} (cos(rln(x)) + 2r*sin(rln(x))) /(1+4r^2)
    where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true).
    """
    first_term = x - 1/2 * log(1 - pow(x, -2)) - log(2*pi)

    second_term = 0
    for r in ls_rz_zero:
        second_term += (cos(r*log(x)) + 2*r*sin(r*log(x))) / (1+4 * pow(r, 2))

    return first_term - 4 * pow(x, 0.5) * second_term


def chebyshev_psi_rz_approximation_2(x, ls_rz_zero):
    """
    Compute the Chebyshev psi function up to x, using
    phi(x) = x - ln(2pi) + 2*sqrt(x)* sum_{r} sin(rln(x)) /r
    where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true).
    """
    first_term = x - log(2*pi)

    second_term = 0
    for r in ls_rz_zero:
        second_term += sin(r*log(x)) / r

    return first_term - 2 * pow(x, 0.5) * second_term


def plot_chebyshev_psi_rz_approximation(x, ls_rz_zero):
    """Plot the Chebyshev psi function up to x."""
    plt.figure(figsize=(10, 6))

    plt.step(range(1, x + 1),
             chebyshev_psi(x)[1:],
             where='post',
             label=r'$\psi(x)$', color='red')

    x_range = np.linspace(2, x, x*20)
    plt.plot(x_range,
             [chebyshev_psi_rz_approximation(i, ls_rz_zero) for i in x_range],
             label=rf'$\psi_0(x)$ with {len(ls_rz_zero)} non-trivial zeros',
             color='blue')

    x_range = np.linspace(2, x, x*20)
    plt.plot(x_range,
             [chebyshev_psi_rz_approximation_2(i, ls_rz_zero) for i in x_range],
             label=rf'$\psi_1(x)$ with {len(ls_rz_zero)} non-trivial zeros',
             color='green')

    # Set x-ticks to be integers
    plt.xticks(np.arange(1, x+1, step=max(x//10, 1)))

    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.title('Chebyshev $\psi$ Prime Counting Function')
    plt.legend()
    # plt.grid(True)

    plt.show()


if __name__ == '__main__':

    # plot_chebyshev_psi(20)

    ls_zeros = []
    with open(r'./imaginary_part_zeros.txt') as f:
        for line in f:
            ls_zeros.append(float(line.strip().split(' ')[1]))

    plot_chebyshev_psi_rz_approximation(20, ls_rz_zero=ls_zeros[:10])
