import unittest
from cmath import exp, pi


class BasicPolynomial:
    def __init__(self, dc_powers):
        # remove 0's
        self.dc_powers = {p: c for p, c in dc_powers.items() if c != 0}

    def __call__(self, x):
        ans = 0
        for p, c in self.dc_powers.items():
            ans += c * (x**p)
        return ans

    def q_eval(self, t, debug=False):
        q = exp(pi*t*complex(0, 1))
        if debug:
            print(q)
        if abs(q.imag) < 1e-14:
            return self.__call__(q.real)
        return self.__call__(q)

    def __add__(self, other):
        dc = self.dc_powers.copy()
        for p, c in other.dc_powers.items():
            dc[p] = dc.get(p, 0) + c
        return BasicPolynomial(dc)

    def __sub__(self, other):
        dc = self.dc_powers.copy()
        for p, c in other.dc_powers.items():
            dc[p] = dc.get(p, 0) - c
        return BasicPolynomial(dc)

    def __mul__(self, other):
        dc = dict()
        for p1, c1 in self.dc_powers.items():
            for p2, c2 in other.dc_powers.items():
                p = p1 + p2
                c = c1*c2
                dc[p] = c + dc.get(p, 0)

        return BasicPolynomial(dc)

    def __eq__(self, other):
        if self.is_constant() and isinstance(other, (int, float)):
            return self.dc_powers.get(1, 0) == other
        return self.dc_powers == other.dc_powers

    def is_constant(self):
        return len(set(self.dc_powers.keys()) - {0}) == 0

    def __str__(self):
        """
        This method returns the string representation of the object. This method is called when print() or str()
        function is invoked on an object.
        This method must return the String object. If we don’t implement __str__() function for a class,
        then built-in object implementation is used that actually calls __repr__() function.
        """
        ret = ''
        for p in sorted(self.dc_powers):
            c = self.dc_powers[p]
            if len(ret) != 0:
                if c > 0:
                    ret += ' + '
                else:
                    ret += ' - '
            if abs(c) == 1:
                if p == 1:
                    ret += 'x'
                elif p == 0:
                    ret += '{}'.format(abs(c))
                elif p < 0:
                    if p == -1:
                        ret += '1/x'
                    else:
                        ret += '1/x^{}'.format(abs(p))
                else:
                    ret += 'x^{}'.format(p)
            else:
                if p == 1:
                    ret += '{}x'.format(abs(c))
                elif p == 0:
                    ret += '{}'.format(abs(c))
                elif p < 0:
                    if p == -1:
                        ret += '{}/x'.format(abs(c))
                    else:
                        ret += '{}/x^{}'.format(abs(c), abs(p))
                else:
                    ret += '{}x^{}'.format(abs(c), p)

        return ret

    def __repr__(self):
        """
        Python __repr__() function returns the object representation.
        It could be any valid python expression such as tuple, dictionary, string etc.
        This method is called when repr() function is invoked on the object, in that case, __repr__()
        function must return a String otherwise error will be thrown.
        """
        return self.dc_powers


class TestPolynomial(unittest.TestCase):

    def test_custom_polynomial(self):
        poly_x = BasicPolynomial({1: 1})  # f(x) = x
        poly_2 = BasicPolynomial({0: 2, 1: -1, 3: 1})  # f(x) = 2 - x + x^3
        poly_3 = BasicPolynomial({-1: -3, 0: 2, 1: -1})  # f(x) = -3/x + 2 - x
        with self.subTest('call x'):
            for i in range(10):
                self.assertEqual(poly_x(i), i)
        with self.subTest('addition'):
            self.assertEqual(BasicPolynomial({0: 2, 3: 1}), poly_x + poly_2)
        with self.subTest('subtraction'):
            self.assertEqual(poly_3 - poly_3, 0)
        with self.subTest('call x for poly 3'):
            f = lambda x: -3/x + 2 - x
            for i in range(1, 10):
                self.assertAlmostEqual(poly_3(i), f(i))
        with self.subTest('string'):
            self.assertEqual(str(poly_3), '3/x + 2 - x')
        with self.subTest('simple multiplication'):
            # x*(2 - x + x^3) = 2x - x^2 + x^4
            self.assertEqual(BasicPolynomial({1: 2, 2: -1, 4: 1}), poly_x*poly_2)
        with self.subTest('multiplication'):
            # (-3/x + 2 - x)*(2 - x + x^3) = -6/x + 7 -4x - 2x^2 +2x^3 - x^4
            self.assertEqual(BasicPolynomial({-1: -6, 0: 7, 1: -4, 2: -2, 3: 2, 4: -1}), poly_3*poly_2)
        with self.subTest('q_eval t=0'):  # t=0 -> q=1
            # -3/x + 2 - x -> -3 + 2 -1 = -2
            self.assertEqual(-2, poly_3.q_eval(0, debug=True))
            # 2 - x + x^3 -> 2 - 1 + 1 = 2
            self.assertEqual(2, poly_2.q_eval(0, debug=True))

        with self.subTest('q_eval t=1'):  # t=1 -> q=-1
            # -3/x + 2 - x -> 3 + 2 +1 = 6
            self.assertEqual(6, poly_3.q_eval(1, debug=True))
            # 2 - x + x^3 -> 2 + 1 - 1 = 2
            self.assertEqual(2, poly_2.q_eval(1, debug=True))
