"""Prime factorization"""
import unittest
from primesieve import primes


def primes_of_n(n, ls_prime=None):
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


class TestPrimes(unittest.TestCase):
    def cases(self, function):
        self.assertEqual(function(5), {5: 1})
        self.assertEqual(function(10), {5: 1, 2: 1})
        self.assertNotEqual(function(4), {5: 1})
        self.assertEqual(function(4), {2: 2})
        self.assertEqual(function(2), {2: 1})
        self.assertEqual(function(32), {2: 5})
        self.assertEqual(function(96), {2: 5, 3: 1})

    def test_no_prime_given(self):
        self.cases(function=primes_of_n)

    def test_prime_given(self):
        self.cases(function=lambda n: primes_of_n(n, ls_prime=primes(100)))
