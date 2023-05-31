"""
See description here: https://www.erdosproblems.com/6

Let d_{n} = p_{n+1} − p_{n}

Are there infinitely many n such that d_{n} < d_{n+1} < d_{n+2} ?

increasing percentage ~= 0.148609 * number of primes
decreasing percentage ~= 0.147695 * number of primes
"""

from primesieve import primes
import matplotlib.pyplot as plt
from math import log


def main():
    n = 10000000
    ls_p = primes(n)
    ls_d = [p2-p1 for p1, p2 in zip(ls_p, ls_p[1:])]
    # print(ls_p)
    # print(ls_d)
    # plt.plot(ls_p, label='prime')
    # plt.plot(ls_d, label='difference')
    # plt.show()

    increasing = 0
    decreasing = 0
    k = 0
    ls_percent_increase = []
    ls_percent_decrease = []
    for d1, d2, d3 in zip(ls_d, ls_d[1:], ls_d[2:]):
        k += 1
        if d1 < d2 < d3:
            print(f'{k=}, increasing {d1=}, {d2=}, {d3=}')
            increasing += 1
        if d1 > d2 > d3:
            print(f'{k=}, decreasing {d1=}, {d2=}, {d3=}')
            decreasing += 1

        ls_percent_increase.append(increasing/k)
        ls_percent_decrease.append(decreasing/k)

    print(f'{n=} has {increasing=}, {decreasing=}')
    plt.plot(ls_percent_increase, label='increasing', color='g')
    plt.plot(ls_percent_decrease, label='decreasing', color='r')
    plt.title('Sequence of differences of primes')
    plt.xlabel('Prime number')
    plt.ylabel('Percentage')
    plt.gca().set_yticklabels(['{:.0f}%'.format(val * 100) for val in plt.gca().get_yticks()])
    plt.legend(loc='lower right')
    plt.show()

    ls_diff_increase = [x1 - x2 for x2, x1 in zip(ls_percent_increase[1:], ls_percent_increase)]
    plt.plot(ls_diff_increase)
    plt.yscale('log')
    plt.xlabel('Prime number')
    plt.ylabel('Log Scale Difference')
    plt.title('Difference of increasing sequence')
    plt.show()


if __name__ == '__main__':
    main()
