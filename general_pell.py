"""Solve general Pell's equation x^2 - D*y^2 = n"""
import unittest
from functools import total_ordering


def is_int(n): return abs(n - int(n)) < 1e-13


@total_ordering
class FieldExtensionD:
    """Numbers of the form a + b*sqrt(D) for a,b integers"""
    def __init__(self, x: int, y: int, d: int):
        self.x = x
        self.y = y
        self.d = d
        self.u_float = x + y * (d ** 0.5)

    def conjugate(self):
        return FieldExtensionD(x=self.x, y=-self.y, d=self.d)

    def norm(self):
        return self*self.conjugate()

    def __mul__(self, other):
        assert self.d == other.d
        f = self.x * other.x + self.y * self.d * other.y
        s = self.y * other.x + self.x * other.y
        return FieldExtensionD(x=f, y=s, d=self.d)

    def __pow__(self, power: int):
        # power is always greater than 0
        if power <= 0:
            raise NotImplementedError('Negative powers are not implemented.')
        result = self
        num_to_mult = power - 1
        while num_to_mult > 0:
            result *= self
            num_to_mult -= 1
        return result

    def get_solutions(self):
        result = self
        while True:
            yield result
            result *= self

    def __str__(self):
        if self.y == 0:
            return f'{self.x}'
        if abs(self.y) == 1:
            if self.y > 0:
                return f'{self.x} + sqrt({self.d})'
            else:
                return f'{self.x} - sqrt({self.d})'
        if self.y < 0:
            return f'{self.x} - {abs(self.y)}*sqrt({self.d})'
        return f'{self.x} + {self.y}*sqrt({self.d})'

    def __repr__(self):
        return str(self)

    def print_generator(self):
        soln = str(self)
        return f"({soln})^k"

    def __eq__(self, other):
        if isinstance(other, FieldExtensionD):
            return (self.x == other.x) and (self.y == other.y) and (self.d == other.d)
        elif isinstance(other, int):
            return (self.x == other) and (self.y == 0)
        else:
            raise NotImplementedError

    def __hash__(self):
        return hash((self.x, self.y, self.d))

    def __lt__(self, other):
        """Returns True if self < other """
        if abs(self.x) < abs(other.x):
            return True
        elif abs(self.x) > abs(other.x):
            return False
        else:
            if self.x > other.x:  # Assume positive x is smaller
                return True
            # only reach here if x's are the same
            if abs(self.y) < abs(other.y):
                return True
            else:
                if self.y > other.y:  # Assume positive y is smaller
                    return True
                else:
                    return False

    @staticmethod
    def remove_duplicates(fundamental_solution, pell_solutions):
        """Filter out duplicate solutions which are a multiple of each other"""
        unique_solutions = pell_solutions.copy()
        removed_solutions = set()  # todo structure this better
        # filter out replicates
        for sol in pell_solutions:
            if sol in removed_solutions:
                continue
            multiple = sol * fundamental_solution
            if multiple in pell_solutions:
                if sol < multiple:
                    unique_solutions.remove(sol)
                    removed_solutions.add(sol)
                else:
                    unique_solutions.remove(multiple)
                    removed_solutions.add(multiple)
        return unique_solutions


class PellEquation:
    """Class to solve Pell Equation x^2 - d*y^2 = 1"""

    def __init__(self, d: int):
        self.d = d
        self.solution = None

    @staticmethod
    def generate_fundamental_solution(d=5):  # todo replace with with convergence continued fraction of sqrt(d)
        """
        Get the fundamental solution to x^2 - d*y^2 = 1, with the smallest positive value of y.
        Args:
            d: <int>, must be non-square
        Returns: <tuple> (x,y) of solution x^2 - d*y^2 = 1, with the smallest positive value of y.
        """
        if is_int(d ** 0.5):
            return None
        y = 1
        while True:
            x2 = 1 + d * y * y
            if is_int(x2 ** 0.5):
                x = int(x2 ** 0.5)
                # only keep the positive values of x and y
                return FieldExtensionD(x=x, y=y, d=d)
            y += 1

    def solve(self):
        if not self.solution:
            self.solution = self.generate_fundamental_solution(d=self.d)
        return self.solution


class GeneralPell(PellEquation):
    """Class to solve Pell Equation x^2 - d*y^2 = n"""

    def __init__(self, d, n):
        super().__init__(d)
        self.n = n
        self.base_solution = super().solve()  # FieldExtensionD

    def solve(self, positive_only=False):
        """
        Get all possibly unique primitive generators of x^2 - d*y^2 = n
        Args:
            positive_only: <bool> specify if positive only solutions should be kept
        Returns: set of FieldExtensionD values of the form (x + y*sqrt(d))
        """
        # need to check |y| <= sqrt(n*u/d)
        u = self.base_solution.u_float
        n, d = self.n, self.d

        abs_y_threshold = int((abs(n) * u / d) ** 0.5)  # y^2 <= n*float_u / d
        if n > 0:
            abs_y_min = 0
        else:
            abs_y_min = int((-n / d) ** 0.5) + 1  # x^2 > 0 => y > sqrt(-n/d)

        ls_tup = set()
        for y in range(abs_y_min, abs_y_threshold + 1):
            x2 = n + d * y * y
            x = int(x2 ** 0.5)
            candidate_solution = FieldExtensionD(x, y, d)
            if candidate_solution.norm() == n:
                if positive_only:
                    # only keep the positive values of x + y*sqrt(d)
                    ls_tup.add(FieldExtensionD(x, y, d))
                    neg_x = FieldExtensionD(-x, y, d)
                    if neg_x.u_float > 0:
                        ls_tup.add(FieldExtensionD(-x, y, d))
                    else:
                        ls_tup.add(FieldExtensionD(x, -y, d))
                    # FieldExtensionD(-x, -y)  # will always be negative
                else:
                    ls_tup.add(FieldExtensionD(x, y, d))
                    ls_tup.add(FieldExtensionD(-x, y, d))
                    ls_tup.add(FieldExtensionD(x, -y, d))
                    ls_tup.add(FieldExtensionD(-x, -y, d))

        return FieldExtensionD.remove_duplicates(self.base_solution, ls_tup)


class TestPell(unittest.TestCase):
    def test_pell_equation(self):
        sol = PellEquation(d=5).solve()
        self.assertEqual((sol.x, sol.y), (9, 4))

    def test_general_pell_equation(self):
        print(GeneralPell(d=6, n=3).solve(positive_only=True))
        print(GeneralPell(d=6, n=-3).solve(positive_only=True))

        self.assertEqual(len(GeneralPell(d=6, n=3).solve(positive_only=True)), 1)
        self.assertEqual(len(GeneralPell(d=6, n=3).solve()), 2)

        print(GeneralPell(d=19, n=36).solve(positive_only=True))
        print(GeneralPell(d=19, n=36).solve(positive_only=False))

    def test_no_general_pell_equation(self):
        """x^2 - 37*y^2 = 11 has no solution"""
        self.assertEqual(len(GeneralPell(d=37, n=11).solve()), 0)

    # # todo confirm solution for n=-1 case
    # def test_neg_general_pell_equation(self):
    #     """x^2 - 5*y^2 = -1"""
    #     self.assertEqual(len(GeneralPell(d=5, n=-1).solve(positive_only=True)), 1)
    #     self.assertEqual(len(GeneralPell(d=5, n=-1).solve(positive_only=False)), 3)

    def test_all_pell_equation(self):
        for d in range(2, 10):
            if is_int(d ** 0.5):
                continue
            with self.subTest('D = {}'.format(d)):
                sol = PellEquation(d=d).solve()
                self.assertEqual(pow(sol.x, 2) - d * pow(sol.y, 2), 1)
