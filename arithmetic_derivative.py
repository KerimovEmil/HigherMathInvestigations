import matplotlib.pyplot as plt


def primes_of_n(n):
    """
    Given an integer n, return a dictionary of prime factors with the keys being the prime number, and the values
    being the multiplicity of that factor.
    """
    if n < 0:
        n = abs(n)
    factors = {}
    nn = n
    i = 2
    while i * i <= nn:
        while nn % i == 0:
            factors[i] = factors.get(i, 0) + 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = 1
    return factors


def log_arithmetic_derivative(n):
    """Computes the log-arithmetic derivative of a number. Currently only works for integers. """
    if n in (0, 1, -1):
        return 0.0
    dc_factorization = primes_of_n(n)
    return sum(m / p for p, m in dc_factorization.items())


def arithmetic_derivative(n):
    """
    Computes the arithmetic derivative of an integer n using exact integer arithmetic:
    n' = sum_{p|n} m * (n // p) where n = prod p^m.
    Rules: 0' = 0, 1' = 0, (-1)' = 0, (-n)' = -n', p' = 1 for prime p.
    """
    if n in (0, 1, -1):
        return 0
    sign = 1 if n > 0 else -1
    abs_n = abs(n)
    dc_factorization = primes_of_n(abs_n)
    return sign * sum(m * (abs_n // p) for p, m in dc_factorization.items())


if __name__ == '__main__':
    assert arithmetic_derivative(5 * 11 * 11 * 11) == 3146
    assert arithmetic_derivative(0) == 0
    assert arithmetic_derivative(1) == 0
    assert arithmetic_derivative(7) == 1
    assert arithmetic_derivative(12) == 16  # 12' = (2^2 * 3)' = 2 * 2 * (12//2) + 1 * (12//3) = 24 + 4 = ? wait: (4*3)' = 4'*3 + 4*3' = 4*3 + 4*1 = 16

    max_n = 1000
    plt.plot(range(max_n), [arithmetic_derivative(i) for i in range(max_n)], 'ro', markersize=2)
    plt.ylabel("Arithmetic Derivative")
    plt.xlabel("Number")
    plt.title("Arithmetic Derivative of Integers")
    plt.show()
