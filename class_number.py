"""
Computing the class number which is the sum of reduced form of a given binary quadratic

external links and references:
    http://oeis.org/A000924
    https://mathworld.wolfram.com/ClassNumber.html
    https://www.springer.com/gp/book/9780387970370
    https://arxiv.org/pdf/1704.00902.pdf
    http://zakuski.utsa.edu/~jagy/indefinite_binary_Buell.pdf
    http://www.numbertheory.org/classnos/
    http://www.numbertheory.org/php/php.html
    http://matwbn.icm.edu.pl/ksiazki/aa/aa83/aa8341.pdf
    odd small class numbers:
      S. Arno, M.L. Robinson, F.S. Wheeler, Imaginary quadratic fields with small odd class number, A
      cta Arith. 83 (1998) 295-330
    even small class numbers:
      From P. Ribenboim, Classical Theory of Algebraic Numbers, p. 636, Springer 2001
      These are squarefree d, not field discriminants, i.e. -d//4 not -d, unless -d mod 4 == 1

"""

from collections import namedtuple
from math import gcd

from numpy import array

# define QuadForm = ax^2 + bxy + cy^2
QuadForm = namedtuple('QuadForm', ['a', 'b', 'c'])


def is_int(n): return abs(n - int(n)) < 1e-13

# reduced form |b| <= a <= c
# discriminant = b^2 - 4ac
# if discriminant = -D
# then only a finite amount of reduced forms. This is the class number
# Proof:
# 4b^2 <= 4ac = b^2 + D
# => 3b^2 <= D
# b <= sqrt(D/3)


class QuadraticForm:

    S = array([[0, -1], [1, 0]])
    T = array([[1, 1], [0, 1]])
    # S = (0, -1)   T = (1, 1)    , since SL2(Z) = <S, T>
    #     (1,  0)       (0, 1)

    def __init__(self, a, b, c):
        assert is_int(a)
        assert is_int(b)
        assert is_int(c)
        self.f = QuadForm(a=a, b=b, c=c)

    def discriminant(self):
        """Return the discriminant"""
        return self.f.b ** 2 - 4 * self.f.a * self.f.c

    def matrix_form(self):
        """
        Return matrix
        A = [ a , b/2]
            [b/2,  c ]
        """
        return array([[self.f.a, self.f.b / 2], [self.f.b / 2, self.f.c]])

    def is_reduced(self):
        """
        Return boolean if quadratic form is in reduced state.
        Reduced is defined as:
        D = discriminant
        if D > 0:
           0 < b < sqrt(D)
           and
           sqrt(D) - b < 2 |a| < sqrt(D) + b
        if D < 0:
           |b| <= a <= c

        # todo see if these are equivalent to
        # -|a| < b <= |a| or 0 <= b <= |a| = |c|

        another definition:
        if |b| ≤ a ≤ c, and b ≥ 0 if either a = c or |b| = a
        """

        # if not self.is_proper():
        #     return False

        def third():
            condition_1 = - abs(self.f.a) < self.f.b <= abs(self.f.a)
            condition_2 = 0 <= self.f.b <= abs(self.f.a)
            condition_3 = abs(self.f.a) == abs(self.f.c)

            return condition_1 or (condition_2 and condition_3)

        def first():
            """
            if D < 0 then |b| <= a <= c
            if D > 0 then 0 < b < sqrt(D) and sqrt(D) - b < 2|a| < sqrt(D) + b
            """
            D = self.discriminant()
            if D < 0:
                return abs(self.f.b) <= self.f.a <= self.f.c
            elif D > 0:
                sr_D = D ** 0.5
                cond_1 = 0 < self.f.b < sr_D
                cond_2 = sr_D - self.f.b < 2 * abs(self.f.a) < sr_D + self.f.b
                return cond_1 and cond_2
            else:
                raise NotImplementedError("D cannot be 0, yet.")

        return third()
        # return first()  # this also works

    def reduce(self):
        """
        Return a reduced form
        """
        pass

    def adjacent(self):
        pass

    def form(self):
        """Returns if form is indefinite or definite"""
        pass

    def __str__(self):
        return f'a: {self.f.a}, b: {self.f.b}, c: {self.f.c}'

    def multiply_by(self, by='T'):
        """
        Returns V'*A*V for V either T or S, as defined in this class
        Args:
            by: <str> 'T' or 'S'

        Returns: an equivalent quadratic form

        """
        if by == 'T':
            v = QuadraticForm.T
        elif by == 'S':
            v = QuadraticForm.S
        else:
            raise NotImplementedError(f'Multiply by {by} is not supported yet')

        new_matrix = v.transpose() @ self.matrix_form() @ v
        assert new_matrix[0][1] == new_matrix[1][0]
        return QuadraticForm(a=int(new_matrix[0][0]), b=int(2 * new_matrix[1][0]), c=int(new_matrix[1][1]))

    def __eq__(self, other):
        if self.f.a == other.f.a and self.f.b == other.f.b and self.f.c == other.f.c:
            return True
        return False

    def is_proper(self):
        """Returns boolean if form is 'proper'. This means not all values are scaled by the same integer"""
        return gcd(gcd(self.f.a, self.f.b), self.f.c) == 1


def is_squarefree(n):
    """
    Check if n is a square-free number, i.e. is divisible by no other perfect square than 1.
    Args:
        n: positive integer to check
    Returns: <boolean> if n is square-free
    """
    for i in range(2, round(n ** 0.5) + 1):
        if n % (i ** 2) == 0:
            return False
    return True


def square_free_sieve(limit):
    """Generator that yields all square free numbers less than limit"""
    a = [True] * limit
    # Needed so we don't mark off multiples of 1^2
    yield 1
    a[0] = a[1] = False
    for i, is_square_free in enumerate(a):
        if is_square_free:
            yield i
            i2 = i * i
            for n in range(i2, limit, i2):
                a[n] = False


def factors_of_n(n):
    """Given an <int> n, return list of tuples of positive integers that multiply to n"""
    assert isinstance(n, int)
    assert n > -1
    ls_out = [(1, n)]
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            ls_out.append((i, n // i))
    return ls_out


def get_reduced_forms(D, debug=False):
    """Given an <int> D, computes all reduced forms of discriminant = -D of the binary quadratic form"""
    ls_potential_reduced_quad_forms = []

    # this filtering was moved out of this function. Not sure which makes more sense. Need to read more theory.
    # cond_1 = ((-D) % 4 == 1) and is_squarefree(D)
    # cond_2 = ((-D) % 4 == 0) and ((-D // 4) % 4 in [2, 3]) and is_squarefree(abs(-D // 4))
    #
    # if not cond_1 and not cond_2:
    #     return []

    abs_b_max = int((D / 3) ** 0.5)
    if debug:
        print(f'max_abs_b: {abs_b_max}')
    for b in range(abs_b_max + 1):
        # temp = 4ac = b^2 + D
        b2 = b ** 2
        temp = b2 + D
        if temp % 4 != 0:
            continue
        else:
            # temp = ac
            temp = temp // 4
            # ac >= b^2
            ls_fac = factors_of_n(temp)
            for fac_tup in ls_fac:
                quad_form = QuadraticForm(a=fac_tup[0], b=b, c=fac_tup[1])
                if quad_form.is_reduced():
                    ls_potential_reduced_quad_forms.append(quad_form)
                    if debug:
                        print(quad_form)
                if b != 0:
                    quad_form = QuadraticForm(a=fac_tup[0], b=-b, c=fac_tup[1])
                    if quad_form.is_reduced():
                        ls_potential_reduced_quad_forms.append(quad_form)
                        if debug:
                            print(quad_form)

    return ls_potential_reduced_quad_forms


def get_class_number(D, debug=False):
    """Given an <int> D, computes the class number of discriminant = -D of the binary quadratic form"""
    ls_reduced = get_reduced_forms(D, debug=debug)

    # sample filtering, which needs to be more advanced
    ls_loop = ls_reduced.copy()
    for potential in ls_loop:
        mult_t = potential.multiply_by(by='T')
        if mult_t in ls_reduced:
            if mult_t != potential:
                ls_reduced.remove(potential)
        mult_s = potential.multiply_by(by='S')
        if mult_s in ls_reduced:
            if mult_s != potential:
                ls_reduced.remove(potential)

    return len(ls_reduced)


def get_negative_class_type(max_n, cn, square_free=False):
    """

    Args:
        max_n:
        cn:
        square_free:

    Returns: sorted <list> of negative discriminants of class number cn up to -max_n

    References:
        https://en.wikipedia.org/wiki/Fundamental_discriminant
    """
    ls_out = []
    print(f'---CLASS NUMBER {cn}---')
    if square_free:
        iter_loop = square_free_sieve(max_n)
    else:
        iter_loop = range(1, max_n)

    for D in iter_loop:

        # The number d denotes a negative fundamental discriminant or, equivalently, the discriminant of an imaginary
        # quadratic number field. In other words, we have either
        # d = 1 (mod 4) and d is square-free,
        # or
        # 4|d and -d/4 = 2 or 3 (mod 4) and d/4 is square free.

        cond_1 = ((-D) % 4 == 1) and is_squarefree(D)
        cond_2 = ((-D) % 4 == 0) and ((-D // 4) % 4 in [2, 3]) and is_squarefree(abs(-D // 4))

        if cond_1 or cond_2:
            if get_class_number(D) == cn:
                if cn % 2 == 0:  # even number
                    if cond_2:
                        fundamental_discriminant = D // 4
                    else:
                        fundamental_discriminant = D
                else:
                    fundamental_discriminant = D
                ls_out.append(fundamental_discriminant)
    ls_out.sort()
    return ls_out


def test_factors_n():
    assert factors_of_n(1) == [(1, 1)]
    assert factors_of_n(2) == [(1, 2)]
    assert factors_of_n(4) == [(1, 4), (2, 2)]


if __name__ == '__main__':
    test_factors_n()

    class_1 = get_negative_class_type(max_n=300, cn=1)
    print(class_1)
    assert class_1 == [3, 4, 7, 8, 11, 19, 43, 67, 163]

    class_2 = get_negative_class_type(max_n=430, cn=2)
    print(class_2)
    assert class_2 == [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    class_3 = get_negative_class_type(max_n=1000, cn=3)
    print(class_3)
    assert class_3 == [23, 31, 59, 83, 107, 139, 211, 283, 307, 331, 379, 499, 547, 643, 883, 907]

    class_4 = get_negative_class_type(max_n=1600, cn=4)
    print(class_4)
    assert class_4 == [14, 17, 21, 30, 33, 34, 39, 42, 46, 55,
                       57, 70, 73, 78, 82, 85, 93, 97, 102, 130,
                       133, 142, 155, 177, 190, 193, 195, 203, 219, 253,
                       259, 291, 323, 355, 435, 483, 555, 595, 627, 667,
                       715, 723, 763, 795, 955, 1003, 1027, 1227, 1243, 1387, 1411, 1435, 1507, 1555]

    assert get_class_number(47) == 5
    assert get_class_number(187) == 2

    # print('---DEBUG 47---')
    # print(get_class_number(47, debug=True))
    # D = 47
    #  a   b   c   works
    #  1   1   12    Y
    #  1  -1   12    X
    #  2   1    6    Y
    #  2  -1    6    Y
    #  3   1    4    Y
    #  3  -1    4    Y

    # note that the solution (1, -1, 12) ~ (1, 1, 12) hence it should not be part of the class number
    # using T = np.array([[1, 1], [0, 1]])
    # we have T' * (1, -1, 12) * T = (1, 1, 12)
    # since
    # (1 0) * (  1   -1/2) * (1 1) = ( 1  1/2)
    # (1 1)   (-1/2   12 )   (0 1)   (1/2  1 )

# import numpy as np
# A = np.array([[1, 1/2], [1/2, 41]])
# S = np.array([[0, -1], [1, 0]])
# T = np.array([[1, 1], [0, 1]])
# A2 = np.array([[1, -1/2], [-1/2, 41]])
# T.transpose() @ A2 @ T


# 187
# x1 = QuadraticForm(a=1, b=1, c=47)
# x2 = QuadraticForm(a=1, b=-1, c=47)
# x3 = QuadraticForm(a=7, b=3, c=7)
# x4 = QuadraticForm(a=7, b=-3, c=7)
# print(x1.multiply_by(by='T'))
# a: 1, b: 3, c: 49
# print(x1.multiply_by(by='S'))
# a: 47, b: -1, c: 1
# print(x2.multiply_by(by='T'))
# a: 1, b: 1, c: 47
# print(x3.multiply_by(by='T'))
# a: 7, b: 17, c: 17
# print(x4.multiply_by(by='T'))
# a: 7, b: 11, c: 11
# print(x3.multiply_by(by='S'))
# a: 7, b: -3, c: 7
# print(x4.multiply_by(by='S'))
# a: 7, b: 3, c: 7
