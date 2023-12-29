"""
Markov Numbers
x,y,z positive integers such that x^2 + y^2 + z^2 = 3xyz

Properties
1. (x,y,z) = (1,1,1) is a solution.
2. If (x,y,z) is a solution then so is (x,z,y),(y,x,z),(y,z,x),(z,x,y),(z,y,x)
3. If (x,y,z) is a solution then so is (x,y,3xy-z)

Proof of 3:
x^2 + y^2 + (3xy-z)^2
= x^2 + y^2 + 9(xy)^2 - 6xyz + z^2
= 3xyz + 9(xy)^2 - 6xyz
= 3xy(z + 3xy - 2z)
= 3xy(3xy - z)
"""
from fractions import Fraction


class MarkovTuple:
    def __init__(self, x: int, y: int, z: int):
        """(x,y,z) such that it satisfies Markov's equation."""
        # self.x = x
        # self.y = y
        # self.z = z
        self.ordered_tup = (min(x, y, z), (x + y + z) - min(x, y, z) - max(x, y, z), max(x, y, z))

        self.unique = max(x, y, z)

    def __str__(self):
        return str(self.ordered_tup)

    def __repr__(self):
        return str(self.ordered_tup)

    def get_next_small_solutions(self):
        x, y, z = self.ordered_tup
        return MarkovTuple(x, z, 3 * x * z - y)

    def get_next_large_solutions(self):
        x, y, z = self.ordered_tup
        return MarkovTuple(y, z, 3 * y * z - x)

    def __eq__(self, other):
        # unproved theorem (uniqueness theorem)
        return self.unique == other.unique

    def __le__(self, other):
        # unproved theorem (uniqueness theorem)
        return self.unique <= other.unique

    def __ge__(self, other):
        # unproved theorem (uniqueness theorem)
        return self.unique >= other.unique

    def get_next_two_solutions(self):
        x, y, z = self.ordered_tup
        a = MarkovTuple(x, z, 3 * x * z - y)
        b = MarkovTuple(y, z, 3 * y * z - x)

        return [a, b]


def get_branching_sol(ls_sol):
    new_ls_sol = ls_sol.copy()
    for sol in ls_sol:

        s1 = sol.get_next_small_solutions()
        s2 = sol.get_next_large_solutions()

        if s1 not in ls_sol:
            new_ls_sol.append(s1)
        if s1 != s2:
            if s2 not in ls_sol:
               new_ls_sol.append(s2)

    return new_ls_sol


def lagrange_sq_numbers(ls_markov):
    """L_n^2 = 9 - 4/m"""
    return [(9 - Fraction(4, m)) for m in ls_markov]


def lagrange_numbers(ls_markov):
    """L_n = sqrt(9 - 4/m)"""
    return [(9 - 4/m)**0.5 for m in ls_markov]


if __name__ == '__main__':
    ls_sol = [MarkovTuple(1, 1, 1)]
    ls_seen_solutions = []
    for _ in range(7):
        new_ls_sol = get_branching_sol(ls_sol=ls_sol)
        solution = new_ls_sol.pop(0)
        print(solution)
        ls_seen_solutions.append(solution)
        ls_sol = new_ls_sol

    print(f'solutions seen: {ls_seen_solutions}')
    print(f'solutions left: {ls_sol}')
    ls_markov = sorted([s.unique for s in (ls_sol + ls_seen_solutions)])
    print(f'Markov Numbers: {ls_markov}')
    print(f'Square LaGrange Numbers: {lagrange_sq_numbers(ls_markov)}')

    # analysis
    # even Markov number mod 32 = 2
    print(f'Even Markov Numbers mod 32: {[m % 32 for m in ls_markov if m%2 == 0]}')
    # odd Markov number mod 4 = 1
    print(f'Odd Markov Numbers mod 4: {[m % 4 for m in ls_markov if m%2 == 1]}')

    print(f'Even Markov Numbers: {[m for m in ls_markov if m % 2 == 0]}')
    # gaps = (0, 4, 12, ...)
    print(f'32n + 2: {[32*n + 2 for n in range(30)]}')
    # index values that are even markov numbers:
    # 0, 1, 6, 19, 202, 1177, 1944, 20188
    # first level of difference,
    # 1, 5, 13, 183, 975, 767, 18244
