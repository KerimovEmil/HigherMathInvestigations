import decimal
from fractions import Fraction
from typing import Callable, Iterator, Union, Optional, List
from math import sin, cos

decimal.getcontext().prec = 100
PI = decimal.Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534'
                     '2117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622'
                     '9489549303819644288109756659334461284756')
E = decimal.Decimal('2.71828182845904523536028747135266249775724709369995')


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

    def fraction(self) -> Fraction:
        """Returns the fraction approximation as defined by self.n_continued_fraction"""
        return self.static_fraction(self.n_continued_fraction)

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
                 min_n=None, ls_an=None, decimal_approx=None):
        self.f = f
        self.f_prime = f_prime
        self.b = 1/100  # can be determined from f'(x) and f''(x) somehow, but 1/100 works for most

        # get initial first few (only two needed) continued fraction values of the number
        if ls_an is not None:
            self.ls_a_n = ls_an
        else:
            min_calculation = min_n if min_n is not None else 2
            self.ls_a_n = ContinuedFraction(decimal_approx).get_continued_fraction(min_calculation)

    @staticmethod
    def get_basic_continued_fraction(x):
        return int(x), 1/(x-int(x))

    def get_ls_a_n(self, max_n: int):

        n = len(self.ls_a_n) + 1
        # get the n-2 convergent
        c_nm2 = cf.static_fraction(self.ls_a_n[:-1])
        # get the n-1 convergent
        c_nm1 = cf.static_fraction(self.ls_a_n)

        while n < max_n:  # get up to max_n terms in continued fraction expansion

            if c_nm2.numerator * c_nm1.denominator - c_nm1.numerator * c_nm2.denominator != (-1) ** n:
                print(f'check failed for n={n}')
                break

            # define the residual
            func_ratio = abs(self.f_prime(c_nm1) / self.f(c_nm1))
            res = func_ratio / (c_nm1.denominator**2) - c_nm2.denominator / c_nm1.denominator

            # setting B >= y_n + 1, ensures that at least one B > y_n
            B = max(self.b * c_nm1.denominator ** 2, c_nm1.denominator + 1)

            m = 0
            while c_nm1.denominator < B:
                m += 1
                a_n, res = self.get_basic_continued_fraction(res)
                self.ls_a_n.append(a_n)
                # set the n-1 and n-2 convergents
                c_nm2 = c_nm1
                c_nm1 = ContinuedFraction.static_fraction(self.ls_a_n)

            n += m
        return self.ls_a_n


def print_continued_fraction(num, n):
    a = ContinuedFraction(num)
    for i in range(n+1):
        cont2 = a.get_continued_fraction(i)
        print(i, cont2, a.fraction(), a.approx(), a.relative_error())
    print('-----------------------------------------------------------------------------------------------')


if __name__ == '__main__':
    # get continued fraction expansion of pi using basic method
    print('get continued fraction of pi using basic method')
    print_continued_fraction(PI, 17)

    # get continued fraction expansion of (e^2-1)/(e^2+1) using basic method
    print('get continued fraction of (e^2-1)/(e^2+1) using basic method')
    cf = ContinuedFraction((E**2 - 1)/(E**2 + 1))
    ls_cf = cf.get_continued_fraction(20)
    print(ls_cf, cf.fraction())

    # define polynomial solution
    print('get continued fraction of 2^(1/3) by defining f(x) = x^3 - 2')
    cf_func_root = ContinuedFractionFunctionRoot(f=lambda x: x**3 - 2, f_prime=lambda x: 3*x**2,
                                                 decimal_approx=2**(1/3))
    print(cf_func_root.get_ls_a_n(max_n=17))

    # # note that this does not work for transcendental numbers like pi, as an error rate needs to be introduced
    # print('get continued fraction of pi by defining f(x) = sin(x)')
    # cf_func_root = ContinuedFractionFunctionRoot(f=sin, f_prime=cos,
    #                                              decimal_approx=3.1415926535, min_n=4)
    # print(cf_func_root.get_ls_a_n(max_n=17))
