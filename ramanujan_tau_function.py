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
    \Delta = 1/1728 * (E_4^3 - E_6^2)
    Args:
        N: highest power to approximate E_4 and E_6
    Returns: BasicPolynomial class
    """
    return 1/1728 * (E(4, N)**3 - E(6, N)**2)


def main():
    approx = 6
    print('E_4 = {}'.format(E(4, approx)))
    print('E_4^3 = {}'.format(E(4, approx)**3))

    print('E_6 = {}'.format(E(6, approx)))
    print('E_6^2 = {}'.format(E(6, approx)**2))

    print('H_10 = {}'.format(E(6, approx)*E(4, approx)))

    delta = discriminant_modular_form(approx)
    print('Delta = {}'.format(delta))


if __name__ == '__main__':
    main()
