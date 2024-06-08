"""
phi(x) = sum_{p^n<=x} log(p) = sum_{k<=x} \Lambda(k)

phi(x) = x - 1/2 ln(1-x^-2) - ln(2pi) + 4*sqrt(x)* sum_{r} (cos(rln(x)) + 2r*sin(rln(x))) /(1+4r^2)
where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true)
"""

import numpy as np
import matplotlib.pyplot as plt
import primesieve
from math import log, sin, cos, pi
from mpmath import li
from utils import memoize, timeit


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


# def mobius_weighted_sum(max_n, ls_f) -> float:
#     partial_sum = 0
#     ls_m = mobius_sieve(max_n)
#
#     for n in range(1, max_n + 1):
#         partial_sum += ls_m[n] / n * ls_f[n]
#
#     return partial_sum


def R(x, ls_primes, max_n=10):
    # return mobius_weighted_sum(max_n=max_n,
    #                            ls_f=[li(pow(x, 1/n), offset=True)
    #                                  for n in range(1, max_n+1)])
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for n in range(1, max_n+1):
        partial_sum += ls_m[n]/n * li(pow(x, 1/n), offset=True)

    return partial_sum


@timeit
def non_trivial_zero_R(x, ls_rz_zero, ls_primes, max_n=10):
    partial_sum = 0
    ls_m = mobius_sieve(max_n + 1, ls_prime=ls_primes)

    for r in ls_rz_zero:
        for n in range(1, max_n + 1):
            partial_sum += 2 * ls_m[n] / n * li(pow(x, complex(0.5, r) / n), offset=True).real

    return partial_sum


@timeit
def prime_counting_rz_approximation(x, ls_rz_zero, ls_primes, max_n=10, max_trivial_zero=10):
    """
    Compute the Chebyshev psi function up to x, using
    pi(x) = R(x) - sum_{p} R(x^p), where p is a zero of the RZ function.
    Assuming the RZ hypothesis
    pi(x) = R(x) - mobius_sum(2*Re(li(x^() sum_{r} (cos(rln(x)) + 2r*sin(rln(x))) /(1+4r^2)
    where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true).
    """
    first_term = R(x, ls_primes=ls_primes, max_n=max_n)
    trivial_zero_term = sum(R(pow(x, -2*m), ls_primes=ls_primes, max_n=max_n)
                            for m in range(1, max_trivial_zero+1))

    non_trivial_term = non_trivial_zero_R(x, ls_rz_zero, ls_primes=ls_primes, max_n=max_n)

    # non_trivial_term = mobius_weighted_sum(max_n=max_n,
    #                                        ls_f=[2*li(pow(x, complex(0.5, r)), offset=True).real
    #                                              for r in ls_rz_zero])
    # non_trivial_term = 0
    # for r in ls_rz_zero:
    #     non_trivial_term += 2 * (cos(r*log(x)) + 2*r*sin(r*log(x)))

    return first_term - trivial_zero_term - non_trivial_term


# def prime_counting_rz_approximation_2(x, ls_rz_zero):
#     """
#     Compute the Chebyshev psi function up to x, using
#     phi(x) = x - ln(2pi) + 2*sqrt(x)* sum_{r} sin(rln(x)) /r
#     where r are the positive imaginary parts of the zeros of the riemann zeta function (if RH is true).
#     """
#     first_term = x - log(2*pi)
#
#     second_term = 0
#     for r in ls_rz_zero:
#         second_term += sin(r*log(x)) / r
#
#     return first_term - 2 * pow(x, 0.5) * second_term


def plot_prime_counting_rz_approximation(max_x, ls_rz_zero, max_trivial_zero=10, max_n=10):
    """
    Plot the prime counting function up to x.
    Args:
        max_x: x-axis for graph
        ls_rz_zero: list of positive imaginary terms of non-trivial RZ zeros
        max_trivial_zero: maximum non-trivial zeros to include
        max_n: how many terms of the mobius inverse summation to use
    """
    plt.figure(figsize=(10, 6))

    plt.step(range(1, max_x + 1),
             [primesieve.count_primes(i) for i in range(1, max_x + 1)],
             where='post',
             label=r'$\pi(x)$', color='red')

    ls_primes = list(primesieve.primes(max_x))
    x_range = np.linspace(2, max_x, max_x*20)
    plt.plot(x_range,
             [prime_counting_rz_approximation(i, ls_rz_zero, ls_primes=ls_primes, max_n=max_n, max_trivial_zero=max_trivial_zero)
              for i in x_range],
             label=rf'$\pi_0(x)$ with {len(ls_rz_zero)} non-trivial zeros and {max_trivial_zero} trivial zeros' +
                   f' and {max_n} terms of mobius inversion',
             color='blue')

    # x_range = np.linspace(2, x, x*20)
    # plt.plot(x_range,
    #          [prime_counting_rz_approximation_2(i, ls_rz_zero) for i in x_range],
    #          label=rf'$\psi_1(x)$ with {len(ls_rz_zero)} non-trivial zeros',
    #          color='green')

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

    # todo fix this currently does not work
    plot_prime_counting_rz_approximation(max_x=100, max_n=50, ls_rz_zero=ls_zeros[:50], max_trivial_zero=0)
