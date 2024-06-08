from math import cos, sin, log


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


def zeta(q: float, N=10) -> complex:
    """returns output of zeta(1/2 + q*i)"""
    return complex(real_zeta(q, N), imag_zeta(q, N))


if __name__ == '__main__':
    print(zeta(14.134725142, N=10))
    print(zeta(14.134725142, N=100))
    print(zeta(14.134725142, N=1000))
    print(zeta(14.134725142, N=10000))
    print(zeta(14.134725142, N=100000))
