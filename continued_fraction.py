import decimal
from fractions import Fraction

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

    def fraction(self) -> Fraction:
        """Returns the fraction approximation as defined by self.n_continued_fraction"""
        ans = Fraction(self.n_continued_fraction[-1])
        for i in self.n_continued_fraction[::-1][1:]:
            ans = i + Fraction(1, ans)
        return ans

    def approx(self) -> float:
        """Returns the approximation as defined by self.n_continued_fraction"""
        return float(self.fraction())

    def relative_error(self) -> float:
        """Returns the relative error of the approximation"""
        return abs(1 - float(self.real_num)/self.approx())


def print_continued_fraction(num, n):
    a = ContinuedFraction(num)
    for i in range(n+1):
        cont2 = a.get_continued_fraction(i)
        print(i, cont2, a.fraction(), a.approx(), a.relative_error())
    print('-----------------------------------------------------------------------------------------------')


if __name__ == '__main__':
    # print_continued_fraction(PI, 10)
    # cf = ContinuedFraction((PI+1)/(PI-1))
    cf = ContinuedFraction((E**2 - 1)/(E**2 + 1))
    ls_cf = cf.get_continued_fraction(10)
    print(ls_cf, cf.fraction())
