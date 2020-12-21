import unittest

# http://oeis.org/A000924
# https://mathworld.wolfram.com/ClassNumber.html
# https://www.springer.com/gp/book/9780387970370
# https://arxiv.org/pdf/1704.00902.pdf
# http://zakuski.utsa.edu/~jagy/indefinite_binary_Buell.pdf


# reduced form |b| <= a <= c
# discriminant = b^2 - 4ac
# if discriminant = -D
# then only a finite amount of reduced forms. This is the class number
# Proof:
# 4b^2 <= 4ac = b^2 + D
# => 3b^2 <= D
# b <= sqrt(D/3)


def factors_of_n(n):
    """Given an <int> n, return list of tuples of positive integers that multiply to n"""
    assert isinstance(n, int)
    assert n > -1
    ls_out = [(1, n)]
    for i in range(2, int(n**0.5)+1):
        if n % i == 0:
            ls_out.append((i, n//i))
    return ls_out


def get_class_number(D, debug=False):
    """Given an <int> D, computes the class number of discriminant = -D of the binary quadratic form"""
    abs_b_max = int((D / 3) ** 0.5)
    if debug:
        print(f'max_abs_b: {abs_b_max}')
    count = 0
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
                if fac_tup[0] >= b and fac_tup[1] >= b:
                    # todo keep list of reduced forms, to check that these new forms are not reduced
                    # example (a,b,c) = (1, -1, 12) ~ (1, 1, 12) through using multiplication of T
                    # T = (1, 1)
                    #     (0, 1)
                    # Where A = ( a , b/2)
                    #         = (b/2,  c )
                    # A' = T^t * A * T
                    count += 1  # todo figure out why this is not +2, since b and -b both work
                    # print(b, b2, temp, fac_tup)
                    if debug:
                        print(f'b: {b}, a: {fac_tup[0]}, c: {fac_tup[1]}')

        # todo implement reduce solutions method which removes multiples of S and T
        # S = (0, -1)   T = (1, 1)    , since SL2(Z) = <S, T>
        #     (1,  0)       (0, 1)

    return count


class TestPrimes(unittest.TestCase):
    def cases(self):
        self.assertEqual(factors_of_n(1), [(1, 1)])
        self.assertEqual(factors_of_n(2), [(1, 2)])
        self.assertEqual(factors_of_n(4), [(1, 4), (2, 2)])


if __name__ == '__main__':
    print('---CLASS NUMBER 1---')
    for i in range(1, 170):
        cn = get_class_number(i)
        if cn == 1:
            print(i)

    print('---CLASS NUMBER 2---')  # todo this list is incorrect
    for i in range(1, 170):
        cn = get_class_number(i)
        if cn == 2:
            print(i)

    print('---DEBUG 47---')
    # this should equal to 5, not 3. Not sure if this is possible to compute.
    print(get_class_number(47, debug=True))
