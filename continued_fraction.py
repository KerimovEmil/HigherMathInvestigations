import unittest
from fractions import Fraction
from typing import Callable, Union


class ContinuedFraction:
    def __init__(self, real_num):
        self.real_num = real_num
        self.n_continued_fraction = None

    def get_continued_fraction(self, n):
        r = self.real_num
        self.n_continued_fraction = []
        for i in range(n+1):
            self.n_continued_fraction.append(int(r))
            r = r - int(r)
            r = 1 / r

        return self.n_continued_fraction

    @staticmethod
    def static_fraction(ls_an) -> Fraction:
        """Returns the fraction approximation as defined by self.n_continued_fraction"""
        ans = Fraction(ls_an[-1])
        for i in ls_an[::-1][1:]:
            ans = i + Fraction(1, ans)
        return ans

    def fraction(self, convergent: int = None) -> Fraction:
        """
        Returns the fraction approximation as defined by self.n_continued_fraction
        Args:
            convergent: if not specified then assuming all of self.n_continued_fraction is used.
        """
        if convergent is None:
            ls_an = self.n_continued_fraction
        else:
            ls_an = self.n_continued_fraction[:convergent]
        return self.static_fraction(ls_an)

    def approx(self) -> float:
        """Returns the approximation as defined by self.n_continued_fraction"""
        return float(self.fraction())

    def relative_error(self) -> float:
        """Returns the relative error of the approximation"""
        return abs(1 - float(self.real_num)/self.approx())


class ContinuedFractionFunctionRoot:
    # implemented algorithm specified by Shiu in
    # https://www.ams.org/journals/mcom/1995-64-211/S0025-5718-1995-1297479-9/S0025-5718-1995-1297479-9.pdf
    # to get continued fraction without requiring a precise decimal expansion
    def __init__(self, f: Callable[[Union[float, Fraction]], float], f_prime: Callable[[Union[float, Fraction]], float],
                 f_prime_ratio: Callable[[Union[float, Fraction]], float] = None,
                 min_n=None, ls_an=None, decimal_approx=None, b=1/100):
        """

        Args:
            f: function such that simple root of f(x)=0 is number we want
            f_prime: derivative of f
            f_prime_ratio: f'/f
            min_n:
            ls_an:
            decimal_approx:
        """
        if f_prime_ratio is None:
            self.f = f
            self.f_prime = f_prime
            self.f_prime_ratio = None
        else:
            self.f_prime_ratio = f_prime_ratio
        self.b = b  # can be determined from f'(x) and f''(x) somehow, but 1/100 works for most

        # get initial first few (only two needed) continued fraction values of the number
        if ls_an is not None:
            self.ls_a_n = ls_an
        else:
            min_calculation = min_n if min_n is not None else 2
            self.ls_a_n = ContinuedFraction(decimal_approx).get_continued_fraction(min_calculation)

    def get_ls_a_n(self, max_n: int):

        n = len(self.ls_a_n) + 1
        # get the n-2 convergent
        c_nm2 = ContinuedFraction.static_fraction(self.ls_a_n[:-1])
        # get the n-1 convergent
        c_nm1 = ContinuedFraction.static_fraction(self.ls_a_n)

        # define the ratio of f'/f
        if self.f_prime_ratio is None:
            def f_prime_ratio(x): return self.f_prime(x) / self.f(x)
        else:
            f_prime_ratio = self.f_prime_ratio

        while n < max_n:  # get up to max_n terms in continued fraction expansion

            if c_nm2.numerator * c_nm1.denominator - c_nm1.numerator * c_nm2.denominator != (-1) ** n:
                print(f'check failed for n={n}')
                break

            # define the residual
            func_ratio = (-1)**(n-1) * f_prime_ratio(c_nm1)  # equivalent to abs(fund_ratio)
            # for non-algebraic functions, we need to ensure precision up to yn^(-4) or y_n^(-2.5), where y_n is the
            # denominator of t_n = x_n / y_n

            res = func_ratio / (c_nm1.denominator**2) - Fraction(c_nm2.denominator, c_nm1.denominator)

            # setting B >= y_n + 1, ensures that at least one B > y_n
            B = max(self.b * c_nm1.denominator ** 2, c_nm1.denominator + 1)
            # B = c_nm1.denominator ** 1.9

            m = 0
            while c_nm1.denominator < B:
                m += 1
                a_n = int(res)
                res = 1 / (res - a_n)
                self.ls_a_n.append(a_n)
                # set the n-1 and n-2 convergents
                c_nm2 = c_nm1
                c_nm1 = ContinuedFraction.static_fraction(self.ls_a_n)
                if n + m > max_n:
                    break

            n += m
        return self.ls_a_n


def print_continued_fraction(num, n):
    a = ContinuedFraction(num)
    for i in range(n+1):
        cont2 = a.get_continued_fraction(i)
        print(i, cont2, a.fraction(), a.approx(), a.relative_error())
    print('-----------------------------------------------------------------------------------------------')


class TestCF(unittest.TestCase):
    def test_cubic_irrational(self):
        cf_func_root = ContinuedFractionFunctionRoot(f=lambda x: x ** 3 - 2,
                                                     f_prime=lambda x: 3 * x ** 2,
                                                     decimal_approx=2 ** (1 / 3))
        with self.subTest('2^(1/3), n=10'):
            self.assertEqual(
                cf_func_root.get_ls_a_n(max_n=11)[:10],
                [1, 3, 1, 5, 1, 1, 4, 1, 1, 8]
                             )
        with self.subTest('2^(1/3), n=70'):
            self.assertEqual(
                cf_func_root.get_ls_a_n(max_n=71)[:70],
                [1, 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, 2, 1, 3, 4, 1, 1, 2, 14, 3, 12, 1, 15,
                 3, 1, 4, 534, 1, 1, 5, 1, 1, 121, 1, 2, 2, 4, 10, 3, 2, 2, 41, 1, 1, 1, 3, 7, 2, 2, 9, 4, 1, 3, 7, 6,
                 1, 1, 2, 2, 9, 3]
                             )

    def test_decimal_method(self):
        """Testing continued fraction using decimal method"""
        import decimal
        decimal.getcontext().prec = 100

        PI = decimal.Decimal(
            '3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534'
            '2117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622'
            '9489549303819644288109756659334461284756')

        # testing pi
        cf = ContinuedFraction(PI)
        with self.subTest('pi, n=17, continued fraction'):
            self.assertEqual(
                cf.get_continued_fraction(17),
                [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2]
                             )
        # test convergent fraction
        with self.subTest('pi, n=2, convergent'):
            self.assertEqual(cf.fraction(convergent=2), Fraction(22, 7))
        with self.subTest('pi, n=17, convergent'):
            self.assertEqual(cf.fraction(), Fraction(2549491779, 811528438))  # all 17

        E = decimal.Decimal('2.71828182845904523536028747135266249775724709369995')

        with self.subTest('(e^2-1)/(e^2+1), n=20, continued fraction'):
            cf = ContinuedFraction((E**2 - 1)/(E**2 + 1))
            self.assertEqual(
                cf.get_continued_fraction(20),
                [0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]
                             )
        with self.subTest('(e^2-1)/(e^2+1), n=20, convergent'):
            self.assertEqual(cf.fraction(),
                             Fraction(371079370602386712421365, 487240307321817004499776))

    def test_quadratic_irrational(self):
        cf_func_root = ContinuedFractionFunctionRoot(f=lambda x: x ** 2 - 2,
                                                     f_prime=lambda x: 2 * x,
                                                     decimal_approx=2 ** (1 / 2))
        with self.subTest('2^(1/2), n=10'):
            self.assertEqual(
                cf_func_root.get_ls_a_n(max_n=11)[:10],
                [1, 2, 2, 2, 2, 2, 2, 2, 2, 2]
                             )
        with self.subTest('2^(1/2), n=50'):
            self.assertEqual(
                cf_func_root.get_ls_a_n(max_n=51)[:50],
                [1] + [2]*49
                             )

    def test_transcendental_pi_using_root(self):
        from math import sin, cos
        # from math import tan

        cf_func_root = ContinuedFractionFunctionRoot(f=sin, f_prime=cos,
                                                     # f_prime_ratio=lambda x: 1/tan(x),
                                                     decimal_approx=3.1415926535)
        with self.subTest('pi, n=13'):
            self.assertEqual(
                cf_func_root.get_ls_a_n(max_n=13)[:13],
                [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
                             )

        # # todo see if we can change B to fix this test
        # # this fails due to the fact that functions like sin(x) are not fractions, hence lose precision
        # # this is an issue for all transcendental numbers like pi, as an error rate needs to be introduced
        # with self.subTest('pi, n=15'):
        #     self.assertEqual(
        #         cf_func_root.get_ls_a_n(max_n=16)[:15],
        #         [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1]
        #                      )


if __name__ == '__main__':
    unittest.main()
