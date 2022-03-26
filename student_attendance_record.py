import unittest

# see derivation of solution here:
# https://massivealgorithms.blogspot.com/2017/04/leetcode-552-student-attendance-record.html

# p(n) = a(n) + p(n-1) + l(n-1)
# l(n) = a(n-1) + p(n-1) + a(n-2) + p(n-2)
# a(n) = a(n-1) + a(n-2) + a(n-3)

# p(1) = 1
# l(1) = 1, l(2) = 3
# a(1) = 1, a(2) = 2, a(3) = 4


def p(n: int) -> int:
    if n == 1:
        return 1
    # return a(n) + p(n-1) + l(n-1)
    # return a(n) + p(n-1) + a(n-2) + p(n-2) + a(n-3) + p(n-3)
    # return a(n) + a(n-2) + a(n-3) + p(n-1) + p(n-2) + p(n-3)
    return a(n) + a(n-2) + a(n-3) + p(n-1) + p(n-2) + p(n-3)


# def l(n: int) -> int:
#     if n == 1:
#         return 1
#     if n == 2:
#         return 3
#     return a(n-1) + p(n-1) + a(n-2) + p(n-2)


def a(n: int) -> int:
    if n == 0:
        return 1
    if n == 1:
        return 1
    if n == 2:
        return 2
    return a(n-1) + a(n-2) + a(n-3)


def t(n: int) -> int:
    """See derivation here: https://brilliant.org/wiki/tribonacci-sequence/"""
    q1 = pow(19 + 3 * pow(33, 0.5), 1/3)
    q2 = pow(19 - 3 * pow(33, 0.5), 1/3)

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


def p_analytical(n: int, c1,c2,c3) -> int:
    # todo compute the exact c1, c2, and c3 that are needed
    q1 = pow(19 + 3 * pow(33, 0.5), 1 / 3)
    q2 = pow(19 - 3 * pow(33, 0.5), 1 / 3)

    x1 = (1 + q1 + q2) / 3
    x2 = 1 / 3 - (q1 + q2) / 6 + complex(0, 1) * (q1 - q2) / (2 * pow(3, 0.5))
    x3 = 1 / 3 - (q1 + q2) / 6 - complex(0, 1) * (q1 - q2) / (2 * pow(3, 0.5))

    b1 = -pow(x1, 2) + 4 * x1 - 1
    b2 = -pow(x2, 2) + 4 * x2 - 1
    b3 = -pow(x3, 2) + 4 * x3 - 1

    # ans = pow(x1, n) / b1 + pow(x2, n) / b2 + pow(x3, n) / b3 + c1*n*pow(x1, n) + c2*n*pow(x2, n) + c3*n*pow(x3, n)
    ans = t(n+1) + c1*n*t(n+1) + c2*n*t(n-1) + c3*n*t(n-2)
    return int(round(ans.real, 0))


class RecurrenceTesting(unittest.TestCase):
    def test_a(self):
        for i in range(1, 10):
            with self.subTest(f'a({i})'):
                self.assertEqual(a_analytical(i), a(i))


if __name__ == '__main__':
    unittest.main()
