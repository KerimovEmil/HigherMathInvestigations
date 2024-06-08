# s(0) = 1
# s(1) = 4

# (k+1)^3 s_{k+1} = 2(2k+1)(5k^2 + 5k +2)s_{k} - 8k(7k^2+1)s_{k-1} + 22k(k-1)*(2k-1)s{k-2}

# [1, 4, 28, 268, 3004, 36784, 476476, 6418192, 88986172, 1261473136, ...]

from utils import memoize


@memoize
def s(k):
    if k == 1:
        return 4
    elif k == 0:
        return 1
    elif k == 2:
        return 28
    else:
        k_sq = k**2
        k_c_inv = (1/k**3)

        first = 2*(2*k-1)*(5*k_sq - 5*k + 2)*s(k-1)
        second = 8*(k-1)*(7*k_sq - 14*k + 8)*s(k-2)
        third = 22*(k-1)*(k-2)*(2*k-3)*s(k-3)

        return int(k_c_inv * (first - second + third))


print([s(i) for i in range(10)])
