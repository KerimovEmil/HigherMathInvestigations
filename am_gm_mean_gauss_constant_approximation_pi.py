from typing import List
import math
import statistics
from functools import reduce


def arithmetic_mean(x: List[float]) -> float:
    return sum(x) / len(x)


def geometric_mean(x: List[float]) -> float:
    product = reduce((lambda i, j: i * j), x)
    return product ** (1/len(x))


def median(x: List[float]) -> float:
    return statistics.median(x)


def agm(x: List[float], tol=1e-4) -> List[float]:
    print(x)
    am = arithmetic_mean(x)
    gm = geometric_mean(x)

    new_x = [am, gm]
    # if max(new_x) - min(new_x) > tol:
    if am - gm > tol:
        return agm(new_x, tol)
    else:
        print(f'am:{am}, gm:{gm}, diff:{am-gm}')
        return new_x


def GMDN(x: List[float], tol=1e-4) -> List[float]:
    print(x)
    am = arithmetic_mean(x)
    gm = geometric_mean(x)
    m = median(x)

    GMDN.num_iter += 1

    new_x = [am, gm, m]
    if max(new_x) - min(new_x) > tol:
        return GMDN(new_x, tol)
    else:
        return new_x


def pi_approx(a=1.0, b=1/math.sqrt(2), t=0.25, p=1, iteration_left=10):
    _pi = pow(a + b, 2) / (4 * t)
    print(f'pi: {_pi}')

    def _update_parameters(a, b, t, p):
        a_n = (a + b) / 2
        b_n = math.sqrt(a * b)
        t_n = t - p * pow(a_n - a, 2)
        p_n = 2 * p
        return a_n, b_n, t_n, p_n

    for i in range(iteration_left):
        _pi = pow(a + b, 2) / (4 * t)
        print(f'iteration: {i}, pi: {_pi}')
        a, b, t, p = _update_parameters(a, b, t, p)


if __name__ == '__main__':
    # print(arithmetic_mean([1, 1, 2, 3, 5]))
    # print(geometric_mean([1, 1, 2, 3, 5]))
    # print(median([1, 1, 2, 3, 5]))

    GMDN.num_iter = 0
    print(GMDN([1, 1, 2, 3, 5], tol=1e-8))  # ~2.089
    print(GMDN.num_iter)

    GMDN.num_iter = 0
    print(f'{GMDN([math.sqrt(2), 1], tol=1e-8)}, iterations: {GMDN.num_iter}')

    # gauss's constant 0.834626841674073186281429732799046808993993013490347002449827370103681992709526411869691160351
    print(1/agm([math.sqrt(2), 1], tol=1e-15)[0])  # ~0.8346268416740731

    # approximation of pi
    pi_approx(iteration_left=5)  # 4 iterations is within 1e-16 of the real value of pi
    print(f'real pi: {math.pi}')
