from math import cos, sin, log
import primesieve


def real_zeta(q: float, N=10) -> float:
    """real part of zeta(1/2 + q*i)"""
    partial_sum = 0
    for n in range(1, N+1):
        partial_sum += cos(q*log(n)) / (n**0.5)
    return partial_sum


def imag_zeta(q: float, N=10) -> float:
    """imaginary part of zeta(1/2 + q*i)"""
    partial_sum = 0
    for n in range(1, N+1):
        partial_sum += sin(q*log(n)) / (n**0.5)
    return -partial_sum


def zeta_original(q: float, N=10) -> complex:
    """returns output of zeta(1/2 + q*i)"""
    return complex(real_zeta(q, N), imag_zeta(q, N))


def zeta(q: float, N=10) -> complex:
    """returns output of zeta(1/2 + q*i)"""
    primes = primesieve.primes(N)
    partial_product = 1

    for p in primes:
        q_ln_p = q * log(p)
        inv_sqrt_p = pow(p, -0.5)
        partial_product *= (1 - inv_sqrt_p * complex(cos(q_ln_p), -sin(q_ln_p)))

    return 1 / partial_product


if __name__ == '__main__':
    # print(zeta_original(14.1347251417346937904572519835625, N=10))
    # print(zeta_original(14.1347251417346937904572519835625, N=1000))
    # print(zeta_original(14.1347251417346937904572519835625, N=100000))

    print(zeta(14.1347251417346937904572519835625, N=10))
    print(zeta(14.1347251417346937904572519835625, N=1000))
    print(zeta(14.1347251417346937904572519835625, N=100000))
    print(zeta(14.1347251417346937904572519835625, N=1000000))
    print(zeta(14.1347251417346937904572519835625, N=10000000))
