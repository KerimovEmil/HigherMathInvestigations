from primesieve import primes
LS_PRIMES = primes(4000)


def primes_of_n(n, ls_prime=LS_PRIMES):
    """
    Given an integer n, return the prime factorization.
    Args:
        n: <int> integer
        ls_prime: <list> optional parameter to specify a list of possible primes
    Returns: <dict> of prime factors with the keys being the prime number, and the values
        being the multiplicity of that factor.
    """
    factors = {}

    if ls_prime is None:
        i = 2
        p = 2

        def next_prime(j):
            return j
    else:
        i = 0
        p = ls_prime[i]

        def next_prime(j):
            return ls_prime[j]

    while p * p <= n:
        while n % p == 0:
            if p not in factors:
                factors[p] = 0
            factors[p] += 1
            n //= p
        i += 1
        p = next_prime(i)

    if n > 1:
        factors[n] = 1
    return factors


class Legendre:
    def __init__(self, num, prime):
        self.p = prime
        self.n = num % self.p

        if self.p == 2:
            raise NotImplementedError("prime equals to 2")
        if self.p not in LS_PRIMES:
            if self.p == 1:
                self.v = 1
            raise NotImplementedError("value is not a prime")
        else:
            self.v = self.value()

    def __mul__(self, other):
        return self.v * other.v

    def __pow__(self, power, modulo=None):
        return self.v ** power

    def __eq__(self, other):
        return self.v == other

    def __repr__(self):
        return str(self.v)

    def value(self):
        if self.n == 0:
            return 0
        elif self.n == -1:
            if self.p % 4 == 1:
                return 1
            else:
                return -1
        elif self.n == 2:
            if self.p % 8 in [1, 7]:
                return 1
            else:
                return -1
        elif self.n**0.5 == int(self.n**0.5):
            return 1
        else:
            if self.n in LS_PRIMES:
                if self.n % 4 == 1 or self.p % 4 == 1:
                    return Legendre(self.p, self.n).v
                else:
                    return -Legendre(self.p, self.n).v
            else:
                if self.n < 0:
                    return Legendre(-1, self.p) * Legendre(-self.n, self.p)
                else:
                    # factor n
                    dc_factors = primes_of_n(self.n)
                    val = 1
                    for p, exp in dc_factors.items():
                        if exp % 2 == 1:
                            val *= Legendre(p, self.p) ** exp
                    return val


if __name__ == '__main__':
    assert Legendre(-1, 5) == 1
    assert Legendre(4, 31) == 1
    assert Legendre(2, 7) == 1
    assert Legendre(24, 31) == -1
    assert Legendre(3411, 3457) == -1
    assert Legendre(7, 7) == 0
