"""B_m = sum_{k=0 to m} (m choose k) B_k"""
import fractions
from scipy.special import comb


def binomial_coefficient(n, r):
    # return combin(n, r)
    return comb(n, r, exact=True)


def combin(n, r):
    """A fast way to calculate binomial coefficients by Andrew Dalke (contrib)."""
    if 0 <= r <= n:
        ntok = 1
        rtok = 1
        for t in range(1, min(r, n - r) + 1):
            ntok *= n
            rtok *= t
            n -= 1
        return ntok // rtok  # bit-wise operation
    else:
        return 0


class BernoulliNumber:
    def __init__(self):
        self.cache_values = {0: fractions.Fraction(1, 1)}  # B_0 = 1

    def get(self, n):
        if n in self.cache_values:
            return self.cache_values[n]

        if n % 2 != 0 and n != 1:
            self.cache_values[n] = fractions.Fraction(0, 1)
            return fractions.Fraction(0, 1)

        b_sum = 0
        for k in range(n):
            # b_sum += combin(n + 1, k) * self.get(k)
            choose_k = fractions.Fraction(binomial_coefficient(n+1, k))
            b_sum += choose_k * self.get(k)
        # b_sum *= -1/(n+1)
        b_sum *= fractions.Fraction(-1, n+1)
        self.cache_values[n] = b_sum
        return b_sum


if __name__ == '__main__':
    b = BernoulliNumber()
    assert b.get(0) == 1
    assert b.get(20) == fractions.Fraction(-174611, 330)
    assert b.get(40) == fractions.Fraction(-261082718496449122051, 13530)

    # Running algo from: https://cs.uwaterloo.ca/journals/JIS/VOL3/KANEKO/AT-kaneko.pdf
    from sympy import symbols

    def f(level, x):
        if level == 0:
            return 1 / x
        else:
            return x * (f(level - 1, x) - f(level - 1, x + 1))

    n = symbols('n')

    for i in range(1, 10):
        term = f(i, n).factor().subs({n: n + 1}).factor()
        bernoulli = term.subs({n: 0})
        print(f'i: {i}, bernoulli: {bernoulli}, poly: {term}')

    # i: 1, bernoulli: 1/2, b*n!: 1, poly: 1/(n + 2)
    # i: 2, bernoulli: 1/6, b*n!: 1, poly: (n + 1)/((n + 2)*(n + 3))
    # i: 3, bernoulli: 0, b*n!: 0, poly: n*(n + 1)/((n + 2)*(n + 3)*(n + 4))
    # i: 4, bernoulli: -1/30, b*n!: -4, poly: (n - 4)*(n + 1)**2/((n + 2)*(n + 3)*(n + 4)*(n + 5))
    # i: 5, bernoulli: 0, b*n!: 0, poly: n*(n + 1)*(n**2 - 13*n - 18)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6))
    # i: 6, bernoulli: 1/42, b*n!: 120, poly: (n + 1)**2*(n**3 - 39*n**2 + 38*n + 120)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7))
    # i: 7, bernoulli: 0, b*n!: 0, poly: n*(n + 1)*(n**4 - 94*n**3 + 371*n**2 + 1546*n + 1200)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8))
    # i: 8, bernoulli: -1/30, b*n!: -12096, poly: (n + 1)**2*(n**5 - 214*n**4 + 2855*n**3 + 3670*n**2 - 9336*n - 12096)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9))
    # i: 9, bernoulli: 0, b*n!: 0, poly: n*(n + 1)*(n**6 - 459*n**5 + 12931*n**4 - 1401*n**3 - 185300*n**2 - 370092*n - 211680)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10))
    # i: 10, bernoulli: 5/66, b*n!: 3024000, poly: (n + 1)**2*(n**7 - 961*n**6 + 54295*n**5 - 309835*n**4 - 1429976*n**3 - 278884*n**2 + 3477360*n + 3024000)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10)*(n + 11))
    # i: 11, bernoulli: 0, b*n!: 0, poly: n*(n + 1)*(n**8 - 1972*n**7 + 199162*n**6 - 2534800*n**5 - 8375831*n**4 + 24013652*n**3 + 125058588*n**2 + 174029040*n + 81648000)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10)*(n + 11)*(n + 12))
    # i: 12, bernoulli: -691/2730, b*n!: -1576143360, poly: (n + 1)**2*(n**9 - 4008*n**8 + 701802*n**7 - 18413688*n**6 + 24101385*n**5 + 509012448*n**4 + 974522828*n**3 - 352544832*n**2 - 2319483456*n - 1576143360)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10)*(n + 11)*(n + 12)*(n + 13))
    # i: 13, bernoulli: 0, b*n!: 0, poly: n*(n + 1)*(n**10 - 8089*n**9 + 2341476*n**8 - 104680266*n**7 + 582069909*n**6 + 5357951367*n**5 + 3331937294*n**4 - 46062600644*n**3 - 134725413000*n**2 - 149105102688*n - 60681519360)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10)*(n + 11)*(n + 12)*(n + 13)*(n + 14))
    # i: 14, bernoulli: 7/6, b*n!: 1525620096000, poly: (n + 1)**2*(n**11 - 16267*n**10 + 7638852*n**9 - 566440902*n**8 + 7995008805*n**7 + 21563873877*n**6 - 193572262402*n**5 - 890138426948*n**4 - 1000972990056*n**3 + 910514773440*n**2 + 2670788937600*n + 1525620096000)/((n + 2)*(n + 3)*(n + 4)*(n + 5)*(n + 6)*(n + 7)*(n + 8)*(n + 9)*(n + 10)*(n + 11)*(n + 12)*(n + 13)*(n + 14)*(n + 15))

    #  4,  5    6    7       8       9        10       11         12            13           14
    # -4, -18, 120, 1200, -12096, -211680, 3024000, 81648000, -1576143360, -60681519360, 1525620096000

# Alternative implementation:

# https://wstein.org/edu/fall05/168/projects/kevin_mcgown/bernproj.pdf
# to get the denominator
# get divisors of m
# add one to them and check if they are prime. If prime multiply them.

# example: m = 2:
# divisors = 1, 2.
# d+1 = 2, 3, both prime. D = 2*3 = 6

# example: m = 4:
# divisors = 1, 2, 4.
# d+1 = 2, 3, 5, all prime. D = 2*3*5 = 30

# example: m = 6:
# divisors = 1, 2, 3, 6.
# d+1 = 2, 3, 4, 7
# d+1 and primes = 2, 3, 7. D = 2*3*7 = 42

# example: m = 8:
# divisors = 1, 2, 4, 8.
# d+1 = 2, 3, 5, 9
# d+1 and primes = 2, 3, 5. D = 2*3*5 = 30

# example: m = 10:
# divisors = 1, 2, 5, 10.
# d+1 = 2, 3, 6, 11
# d+1 and primes = 2, 3, 11. D = 2*3*11 = 66

# example: m = 14:
# divisors = 1, 2, 7, 14.
# d+1 = 2, 3, 8, 15
# d+1 and primes = 2, 3. D = 2*3 = 6
