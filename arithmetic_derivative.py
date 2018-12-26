import matplotlib.pyplot as plt


def primes_of_n(n):
    """
    Given an integer n, return a dictionary of prime factors with the keys being the prime number, and the values
    being the multiplicity of that factor.
    """
    factors = {}
    nn = n
    i = 2
    while i * i <= nn:
        while nn % i == 0:
            if i not in factors:
                factors[i] = 0
            factors[i] += 1
            nn //= i
        i += 1
    if nn > 1:
        factors[nn] = 1
    return factors


def log_arithmetic_derivative(n):
    """Computes the log-arithmetic derivative of a number. Currently only works for integers. """
    dc_factorization = primes_of_n(n)
    return sum([m/p for p, m in dc_factorization.items()])


def arithmetic_derivative(n):
    """Computes the arithmetic derivative of a number. Currently only works for integers. """
    return n * log_arithmetic_derivative(n)


if __name__ == '__main__':
    assert arithmetic_derivative(5*11*11*11) == 3146

    x = 100000
    # plt.plot(range(x), [arithmetic_derivative(x) for x in range(x)])
    plt.plot(range(x), [arithmetic_derivative(x) for x in range(x)], 'ro')
    plt.ylabel("Arithmetic Derivative")
    plt.xlabel("Number")
    plt.show()

