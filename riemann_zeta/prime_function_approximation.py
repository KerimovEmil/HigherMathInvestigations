"""
References: http://ism.uqam.ca/~ism/pdf/Hutama-scientific%20report.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
import primesieve
from math import log, sin, cos, pi, atan
from mpmath import li
from utils import memoize, timeit


def f_li(x: complex) -> complex:
    return li(x, offset=True)
    # return x / log(x)


@memoize
def mobius_sieve(n: int, ls_prime: list[int]) -> list[int]:
    """
    Returns a list of all mobius function values.
    mobius(n) = 1 if i is square-free with even number of primes,
               -1 if odd number,
                0 if contains square
    """
    ls_m = [1]*n
    for p in ls_prime:
        ls_m[p:n:p] = [-1 * x for x in ls_m[p:n:p]]
        p2 = p ** 2
        ls_m[p2:n:p2] = [0] * ((n-1)//p2)  # len(ls_m[p2:n:p2]) == (n-1)//p2
    return ls_m


def R(x, ls_primes, max_n=10):
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for n in range(1, max_n+1):
        mu = ls_m[n]
        if mu != 0:
            partial_sum += mu / n * f_li(pow(x, 1/n))

    return partial_sum


@timeit
def non_trivial_zero_R_1(x, ls_rz_zero, ls_primes, max_n=10):
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for r in ls_rz_zero:
        for n in range(1, max_n + 1):
            mu = ls_m[n]
            if mu != 0:
                partial_sum += 2 * mu / n * f_li(pow(x, complex(0.5, r) / n)).real

    return partial_sum


@timeit
def non_trivial_zero_R_2(x, ls_rz_zero, ls_primes, max_n=10):
    """Simplifying using li(x) ~ x/ln(x)"""
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for n in range(1, max_n + 1):
        mu = ls_m[n]
        if mu != 0:
            temp_n_partial = mu * pow(x, 1/(2*n))
            temp_r_partial = 0
            for r in ls_rz_zero:
                arg = log(x) * r / n
                temp_r_partial += (cos(arg) + 2*r*sin(arg)) / (1 + 4 * pow(r, 2))

            partial_sum += temp_n_partial * temp_r_partial

    return partial_sum * 4 / log(x)


@timeit
def non_trivial_zero_R_3(x, ls_rz_zero, ls_primes, max_n=10):
    """Simplifying using li(x) ~ x/ln(x) and large x simplifications"""
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for n in range(1, max_n + 1):
        mu = ls_m[n]
        if mu != 0:
            temp_n_partial = mu * pow(x, 1/(2*n))
            temp_r_partial = 0
            for r in ls_rz_zero:
                arg = log(x) * r / n
                temp_r_partial += sin(arg) / (2 * r)

            partial_sum += temp_n_partial * temp_r_partial

    return partial_sum * 4 / log(x)

@timeit
def trivial_zero_R_1(x, ls_primes, max_n=10, num_trivial_zero=10):
    return sum(R(pow(x, -2*m), ls_primes=ls_primes, max_n=max_n) for m in range(1, num_trivial_zero+1))


@timeit
def trivial_zero_R_2(x, ls_primes, max_n=10):
    """Simplifying using li(x) ~ x/ln(x)"""
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for n in range(1, max_n + 1):
        mu = ls_m[n]
        if mu != 0:
            partial_sum += mu * log(1 - pow(x, -2/n))

    return partial_sum / (2 * log(x))


@timeit
def trivial_zero_R_3(x):
    """From li(x) = Ei(ln(x)) we get 1/ln(x) - arctan(pi/ln(x))/pi"""
    return 1/log(x) - atan(pi/log(x)) / pi


@timeit
def prime_counting_rz_approximation(x, ls_rz_zero, ls_primes, max_n=10, num_trivial_zero=10):
    """
    Compute the Chebyshev psi function up to x, using
    pi(x) = R(x) - sum_{p} R(x^p), where p is a zero of the RZ function.
    Assuming the RZ hypothesis
    pi(x) = R(x) - mobius_sum(2*Re(li(x^()
    where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true).
    """
    first_term = R(x, ls_primes=ls_primes, max_n=max_n)

    # trivial_zero_term = trivial_zero_R_1(x, ls_primes, max_n=max_n, num_trivial_zero=max_trivial_zero)
    # trivial_zero_term = trivial_zero_R_2(x, ls_primes, max_n=max_n)
    trivial_zero_term = trivial_zero_R_3(x)

    # non_trivial_term = non_trivial_zero_R_1(x, ls_rz_zero, ls_primes=ls_primes, max_n=max_n)
    # non_trivial_term = non_trivial_zero_R_2(x, ls_rz_zero, ls_primes=ls_primes, max_n=max_n)
    non_trivial_term = non_trivial_zero_R_3(x, ls_rz_zero, ls_primes=ls_primes, max_n=max_n)

    return first_term - trivial_zero_term - non_trivial_term


def plot_prime_counting_rz_approximation(max_x, ls_rz_zero, num_trivial_zero=10, max_n=None):
    """
    Plot the prime counting function up to x.
    Args:
        max_x: x-axis for graph
        ls_rz_zero: list of positive imaginary terms of non-trivial RZ zeros
        num_trivial_zero: number of non-trivial zeros to include
        max_n: how many terms of the mobius inverse summation to use
    """
    plt.figure(figsize=(10, 6))

    if max_n is None:
        # max_n = int(log(max_x) / log(2)) + 1
        max_n = int(log(max_x) - log(2)) + 2

    plt.step(range(1, max_x + 1),
             [primesieve.count_primes(i) for i in range(1, max_x + 1)],
             where='post',
             label=r'$\pi(x)$', color='red')

    ls_primes = list(primesieve.primes(max_n + 1))
    x_range = np.linspace(2, max_x, max_x*20)
    plt.plot(x_range,
             [prime_counting_rz_approximation(i, ls_rz_zero, ls_primes=ls_primes, max_n=max_n, num_trivial_zero=num_trivial_zero)
              for i in x_range],
             label=rf'$\pi_0(x)$ with {len(ls_rz_zero)} non-trivial zeros and {num_trivial_zero} trivial zeros' +
                   f' and {max_n} terms of mobius inversion',
             color='blue')

    # Set x-ticks to be integers
    plt.xticks(np.arange(1, max_x+1, step=max(max_x//10, 1)))

    plt.xlabel('x')
    plt.ylabel(r'$\pi(x)$')
    plt.title(r'$\pi(x)$ Prime Counting Function')
    plt.legend()

    plt.show()


if __name__ == '__main__':

    ls_zeros = []
    with open(r'./imaginary_part_zeros.txt') as f:
        for line in f:
            ls_zeros.append(float(line.strip().split(' ')[1]))

    plot_prime_counting_rz_approximation(max_x=100, max_n=None, ls_rz_zero=ls_zeros[:100], num_trivial_zero='all')
