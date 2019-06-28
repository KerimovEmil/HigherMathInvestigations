"""Guess's the positive integer coefficient of polynomial f(x), with only two questions"""


def format_equation(ls_coef):
    """
    Given a list of coefficients outputs a string of the polynomial

    Args:
        ls_coef: <list> a list of coefficients (in ascending order)

    Returns: <str> polynomial

    """
    d = len(ls_coef)
    o_str = 'f(x) = '
    for i in range(d - 1, 0, -1):
        if ls_coef[i] != 0:
            o_str += str(ls_coef[i]) + 'x^' + str(i) + ' + '
    o_str += str(ls_coef[0])

    return o_str


def find_int_poly(f_1, f_f_1):
    """
    Given the value of the integer polynomial at x=1 and at x = f(1) + 1, this function outputs the original integer
    polynomial
    Args:
        f_1: f(1)
        f_f_1: f(f(1) + 1)

    Returns: the string of f

    """

    coeff = []
    N = f_1 + 1
    yi = f_f_1

    while yi > 0:
        ai = yi % N
        coeff.append(ai)
        yi = (yi - ai) // N

    return format_equation(coeff)


if __name__ == '__main__':
    def f(x):
        """Arbitrary make function f(x)"""
        return 678 * (x ** 5) + 222 * (x ** 3) + 1 * x ** 2 + 0

    def g(x):
        """Arbitrary make function f(x)"""
        return 6412 * (x ** 7) + 142 * (x ** 4) + 132 * (x ** 3) + 1 * x ** 2 + 5*x

    assert find_int_poly(f(1), f(f(1) + 1)) == "f(x) = 678x^5 + 222x^3 + 1x^2 + 0"
    assert find_int_poly(g(1), g(g(1) + 1)) == "f(x) = 6412x^7 + 142x^4 + 132x^3 + 1x^2 + 5x^1 + 0"

    # Question 1, what is the value of the polynomial at n=1?
    # p = f(1)
    first_answer = 901  # sample answer

    # Question 2, what is the value of the polynomial at n= f(1) + 1?
    # poly = f(N)
    second_answer = 404820555383370676  # sample answer

    polynomial = find_int_poly(first_answer, second_answer)

    print(first_answer, second_answer, polynomial)



