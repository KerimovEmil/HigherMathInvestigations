"""
https://en.wikipedia.org/wiki/Divisor_function
σ_{x}(n), for a real or complex number x, is defined as the sum of the xth powers of the positive divisors of n.
It can be expressed in sigma notation as

σ_{x}(n) = sum_{d|n} d^x

for x !=0:
σ_{x}(n) = prod_{i=1}^{r} (p_i^((a_i+1)*x) - 1) / (p_i^x - 1)

Recurrence relation:
# todo

"""
from prime_factorization import primes_of_n
from functools import reduce
import unittest


def multiply(ls):
    """Multiply each element of the iterator"""
    return reduce(lambda a, b: a * b, ls)


def num_of_divisors(n):
    """
    Return the number of positive divisors of n.
    e.g. sigma(12) = 6
    """
    dc_primes = primes_of_n(n)
    return multiply(v + 1 for v in dc_primes.values())


def divisors(x, n):
    """return the divisor function σ_{x}(n) = sum_{d|n} d^x"""
    if n == 1:
        return 1
    if x == 0:
        return num_of_divisors(n)

    # for x !=0:
    # σ_{x}(n) = prod_{i=1}^{r} (p_i^((a_i+1)*x) - 1) / (p_i^x - 1)
    dc_primes = primes_of_n(n)
    return multiply(
        (p**((e + 1)*x) - 1) / (p**x - 1)
        for p, e in dc_primes.items())


class TestDivisors(unittest.TestCase):

    def test_num_of_divisors(self):
        self.assertEqual(num_of_divisors(12), 6)
        self.assertEqual(divisors(0, 12), 6)

    def test_other_x(self):
        self.assertEqual(divisors(0, 12), 6)

    def test_x_1(self):
        ls_outputs = [1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12]
        for i, n in enumerate(range(1, 12)):
            self.assertEqual(divisors(x=1, n=n), ls_outputs[i])

    def test_x_2(self):
        ls_outputs = [1, 5, 10, 21, 26, 50, 50, 85, 91, 130, 122]
        for i, n in enumerate(range(1, 12)):
            self.assertEqual(divisors(x=2, n=n), ls_outputs[i])

    def test_x_3(self):
        ls_outputs = [1, 9, 28, 73, 126, 252, 344, 585, 757, 1134, 1332]
        for i, n in enumerate(range(1, 12)):
            self.assertEqual(divisors(x=3, n=n), ls_outputs[i])

    def test_x_4(self):
        ls_outputs = [1, 17, 82, 273, 626, 1394, 2402, 4369, 6643, 10642, 14642]
        for i, n in enumerate(range(1, 12)):
            self.assertEqual(divisors(x=4, n=n), ls_outputs[i])
