from polynomial import BasicPolynomial
from bernoulli_numbers import BernoulliNumber
from divisor import divisors


def E(k, N):
    """
    E_k = 1 - 2*k / B_k * sum_{n=1}^{N} sigma_{k-1}(n) q^n
    Args:
        k: weight
        N: Get up to the q^N power approximation

    Returns:

    """
    dc = {0: 1}
    coeff = int(-2*k / BernoulliNumber().get(k))

    for i in range(1, N+1):
        dc[i] = coeff*int(divisors(k-1, i))
    return BasicPolynomial(dc, 'q')


def discriminant_modular_form(N):
    """
    /Delta = 1/1728 * (E_4^3 - E_6^2)
    Args:
        N: highest power to approximate E_4 and E_6
    Returns: BasicPolynomial class
    """
    return 1/1728 * (E(4, N)**3 - E(6, N)**2)


def j_invariant(N):
    """
    j = E_4^3 / Delta
    Args:
        N: highest power to approximate E_4 and E_6
    Returns: BasicPolynomial class
    """
    return E(4, N) ** 3 / discriminant_modular_form(N)


def print_ramanujan_tau_function(delta, N):
    for i in range(1, N+1):
        print('T({}) = {}'.format(i, delta.dc_powers.get(i, 0)))


def main():
    approx = 10
    # print('E_4 = {}'.format(E(4, approx)))
    print('E_4^3 = {}'.format(E(4, approx)**3))
    print('1 / E_4^3 = {}'.format((E(4, approx)**3).invert(approx)))

    print('E_6 = {}'.format(E(6, approx)))
    print('E_6^2 = {}'.format(E(6, approx)**2))

    delta = discriminant_modular_form(approx)
    print('Delta = {}'.format(delta))
    # print_ramanujan_tau_function(delta, approx)

    # j_func = j_invariant(approx)
    # print('j = {}'.format(j_func))


if __name__ == '__main__':
    main()
