import decimal

decimal.getcontext().prec = 100
PI = decimal.Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534'
                     '2117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622'
                     '9489549303819644288109756659334461284756')


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


def print_continued_fraction(num, n):
    a = ContinuedFraction(num)
    for i in range(n+1):
        cont2 = a.get_continued_fraction(i)
        print(i, cont2)
    print('-----------------------------------------------------------------------------------------------')


print_continued_fraction(PI, 100)
