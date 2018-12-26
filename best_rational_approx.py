from fractions import Fraction
import decimal

decimal.getcontext().prec = 100
PI = decimal.Decimal('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534'
                     '2117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622'
                     '9489549303819644288109756659334461284756')


class AdditionalFraction(Fraction):
    def mediant(self, b):
        """mediant(a, b)"""
        return AdditionalFraction(self.numerator + b.numerator, self.denominator + b.denominator)

    def value(self):
        return float(self.numerator/self.denominator)


class ContinuedFraction:
    def __init__(self, real_num):
        self.real_num = real_num
        self.n_continued_fraction = None

    def get_n_approximation(self, n):
        steps = int(self.real_num)
        starting = AdditionalFraction(steps + 1, 1)
        step_size = AdditionalFraction(steps, 1)
        upper = True

        self.n_continued_fraction = [steps]

        for i in range(n):
            starting, step_size, upper, steps = self.get_next_approx(starting, step_size, upper)
            self.n_continued_fraction.append(steps)

        return step_size

    def get_continued_fraction(self):
        return self.n_continued_fraction

    def get_next_approx(self, starting, step_size, upper):
        steps = 0
        while True:
            previous = starting
            starting = starting.mediant(step_size)
            steps += 1
            if upper:
                if starting.value() < self.real_num:
                    break
            else:
                if starting.value() > self.real_num:
                    break
        new_starting = starting
        new_step_size = previous
        upper = not upper
        return new_starting, new_step_size, upper, steps


def print_continued_frac(num, n):
    a = ContinuedFraction(num)
    for i in range(n+1):
        frac = a.get_n_approximation(i)
        cont = a.get_continued_fraction()
        print(i, frac, frac.value(), cont, frac.value() - float(num))


print_continued_frac(PI, 13)
print_continued_frac(decimal.Decimal(2**0.5), 13)
