import unittest
from functools import lru_cache

# see derivation of solution here:
# https://massivealgorithms.blogspot.com/2017/04/leetcode-552-student-attendance-record.html

# p(n) = a(n-1) + p(n-1) + l(n-1)
# l(n) = a(n-1) + p(n-1) + a(n-2) + p(n-2)
# a(n) = a(n-1) + a(n-2) + a(n-3)

# p(1) = 1
# l(1) = 1, l(2) = 3
# a(1) = 1, a(2) = 2, a(3) = 4

# these relations can be summarized as
# a(n) = a(n-1) + a(n-2) + a(n-3), a(1)=1, a(2)=2, a(3)=4
# p(n) = p(n-1) + p(n-2) + p(n-3) + a(n), p(0)=0, p(1)=1, p(2)=4

# note that total(n) = a(n) + p(n) + L(n)
# this simplifies to total(n) = p(n+1)


@lru_cache(maxsize=None)
def p(n: int) -> int:
    if n == 1:
        return 1
    return a(n-1) + p(n-1) + L(n-1)


@lru_cache(maxsize=None)
def L(n: int) -> int:
    if n == 1:
        return 1
    if n == 2:
        return 3
    return a(n-1) + p(n-1) + a(n-2) + p(n-2)


@lru_cache(maxsize=None)
def a(n: int) -> int:
    if n == 0:
        return 1
    if n == 1:
        return 1
    if n == 2:
        return 2
    return a(n-1) + a(n-2) + a(n-3)


def t(n: int) -> int:
    """See derivation here: https://github.com/KerimovEmil/MathLatexDocs/blob/main/pdfs/tribonacci_numbers_main.pdf """
    q1 = pow(19 + 3 * pow(33, 0.5), 1/3)
    q2 = pow(19 - 3 * pow(33, 0.5), 1/3)

    # x1, x2, x3 = alpha, beta, gamma = roots of x^3 = x^2 + x + 1
    x1 = (1 + q1 + q2) / 3
    x2 = 1/3 - (q1 + q2) / 6 + complex(0, 1) * (q1 - q2) / (2 * pow(3, 0.5))
    x3 = 1/3 - (q1 + q2) / 6 - complex(0, 1) * (q1 - q2) / (2 * pow(3, 0.5))

    b1 = -pow(x1, 2) + 4*x1 - 1
    b2 = -pow(x2, 2) + 4*x2 - 1
    b3 = -pow(x3, 2) + 4*x3 - 1

    ans = pow(x1, n) / b1 + pow(x2, n) / b2 + pow(x3, n) / b3

    return int(round(ans.real, 0))


def a_analytical(n: int) -> int:
    """a(n) = t(n+1)"""
    return t(n+1)


def p_analytical(n: int) -> int:
    ans = 1/22 * (10*n * t(n+1) + (7+5*n)*t(n) + (3+3*n)*t(n-1))
    return int(round(ans, 0))


class RecurrenceTesting(unittest.TestCase):
    def test_a(self):
        for i in range(1, 10):
            with self.subTest(f'a({i})'):
                self.assertEqual(a_analytical(i), a(i))

    def test_p(self):
        for i in range(1, 10):
            with self.subTest(f'p({i})'):
                self.assertEqual(p_analytical(i), p(i))


if __name__ == '__main__':
    unittest.main()

# solving for recurrence solution for p(n)
# from sympy import Symbol, Eq, solve
#
# a1 = Symbol('a1')
# a2 = Symbol('a2')
# a3 = Symbol('a3')
# a4 = Symbol('a4')
# b1 = Symbol('b1')
# b2 = Symbol('b2')
# b3 = Symbol('b3')
# b4 = Symbol('b4')
#
#
# def t(n: int) -> int:
#     if n == -2:
#         return 1
#     if n == -1:
#         return 0
#     if n == 0:
#         return 0
#     if n == 1:
#         return 1
#     if n == 2:
#         return 1
#     return t(n-1) + t(n-2) + t(n-3)
#
# def p(n: int) -> int:
#     if n == 1:
#         return 1
#     if n == 2:
#         return 3
#     if n == 3:
#         return 8
#     return t(n) + t(n-1) + t(n-2) + p(n-1) + p(n-2) + p(n-3)

#
# ls_eq = []
# for n in range(16):
#     eq = Eq(p(n), (a1 + b1 * n) * t(n+1) + (a2 + b2 * n) * t(n) + (a3 + b3 * n) * t(n-1) + (a4 + b4 * n) * t(n-2))
#     ls_eq.append(eq)
#
# ans = solve(ls_eq, (a1, a2, a3, a4, b1, b2, b3, b4))
