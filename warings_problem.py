# if 2^k frac( (3/2)^k ) + floor( (3/2)^k ) <= 2^k for all k then
# g(k) = 2^k + floor( (3/2)^k ) - 2

# this eq has been proved valid for k < 471,600,000

# frac((3/2)^k) = (3/2)^k - floor((3/2)^k)

# 2^k ((3/2)^k - floor( (3/2)^k )) + floor( (3/2)^k ) <= 2^k
# 3^k + (1-2^k) * floor( (3/2)^k )) <= 2^k
# 3^k - 2^k <= (2^k - 1) * floor( (3/2)^k ))

# (3/2)^k - 1 <= (1 - 2^-k) * floor( (3/2)^k ))

# floor( (3/2)^k )) => ((3/2)^k - 1) / (1 - 2^-k)
# floor( (3/2)^k )) / (3/2)^k => (1- (2/3)^k) / (1 - 2^-k)
# floor((3/2)^k )) / (3/2)^k => (6^k- 4^k) / (6^k - 3^k)


import matplotlib.pyplot as plt
import math

# m = 100

# # plotting the difference
# # diff = [int((3/2)**i) * (2/3)**i - (6**i - 4**i) / (6**i - 3**i) for i in range(2, m)]
# diff = [int(math.pow(1.5, i)) * math.pow(2/3, i)
#         -
#         (1 - math.pow(2/3, i)) / (1 - math.pow(2, -i))
#         for i in range(2, m)]
# plt.plot(diff)
# plt.yscale('log')
# plt.ylabel('Inequality greater than 0')
# plt.xlabel('k')
# plt.legend(['floor((3/2)^k )) / (3/2)^k - (6^k- 4^k) / (6^k - 3^k)'])

# # plotting the floor((3/2)^k )) / (3/2)^k
# floor_ratio = [1 - int((3/2)**i) / (3/2)**i for i in range(2, m)]
# plt.plot(floor_ratio)
# plt.yscale('log')
# plt.ylabel('Behaviour of floor((3/2)^k )) / (3/2)^k')
# plt.xlabel('k')
# plt.legend(['1 - floor((3/2)^k )) / (3/2)^k '])

# # plotting the (6^k- 4^k) / (6^k - 3^k)
# power_ratio = [1 - (6**i - 4**i) / (6**i - 3**i) for i in range(2, m)]
# plt.plot(power_ratio)
# plt.yscale('log')
# plt.ylabel('Behaviour of (6^k- 4^k) / (6^k - 3^k)')
# plt.xlabel('k')
# plt.legend(['1 - (6^k- 4^k) / (6^k - 3^k)'])

from decimal import *

# floor((3/2)^k )) / (3/2)^k => (6^k- 4^k) / (6^k - 3^k)

# floor((3/2)^k) * (2^k -1) / (3^k - 2^k) - 1 > 0
# floor((3/2)^k) * (2^k -1) / (3^k - 2^k) > 1
# floor((3/2)^k) * (2^k - 1) > (3^k - 2^k)


def diffa(a, b):
    m = 1750
    getcontext().prec = 1000
    a_dec = Decimal(a)
    b_dec = Decimal(b)
    return [int(a_dec**i / b_dec**i) * (b_dec**i - 1) / (a_dec**i - b_dec**i) - 1 for i in range(2, m)]


def ls_diff(up_to=1750, precision=1000):
    """
    Return list of
    floor((3/2)^k) * (2^k -1) / (3^k - 2^k) - 1
    for k from 1 to up_to

    If this is always positive then floor((3/2)^k) * (2^k - 1) > (3^k - 2^k)
    """
    getcontext().prec = precision
    dec_3 = Decimal(3)
    dec_2 = Decimal(2)
    return [int(dec_3**i / dec_2**i) * (dec_2**i - 1) / (dec_3**i - dec_2**i) - 1 for i in range(2, up_to)]


def ls_diff_2(up_to=1750, precision=1000):
    """
    Return list of
    floor((3/2)^k )) / (3/2)^k - (1- (2/3)^k) / (1 - 2^-k)
    for k from 1 to up_to

    If this is always positive then floor((3/2)^k )) / (3/2)^k => (6^k- 4^k) / (6^k - 3^k)
    """
    getcontext().prec = precision
    dec_2 = Decimal(2)
    dec_1_5 = Decimal(1.5)
    return [int(dec_1_5**i) / (dec_1_5**i) - (1 - dec_1_5**(-i)) / (1 - dec_2**(-i)) for i in range(2, up_to)]


def ls_log_diff(up_to=1750, precision=1000):
    """
    Return list of
    log( floor((3/2)^k) ) + log((2^k -1)) - log(3^k - 2^k) > 0
    for k from 1 to up_to

    If this is always positive then floor((3/2)^k) * (2^k - 1) > (3^k - 2^k)
    """
    getcontext().prec = precision
    dec_3 = Decimal(3)
    dec_2 = Decimal(2)
    dec_1_5 = dec_3 / dec_2

    ls_out = []

    # math.log(dec_3 ** k - dec_2 ** k)
    # math.log(dec_2 ** k * (dec_1.4 ** k - 1))
    # k * math.log(dec_2) + math.log((dec_1.5 ** k - 1))

    for k in range(2, up_to):
        three_halfs_raised = dec_1_5**k
        n = math.log(int(three_halfs_raised) * (dec_2**k - 1)/(three_halfs_raised - 1)) - k * math.log(dec_2)
        ls_out.append(n)
    return ls_out


diff = ls_diff(1800, precision=1000)
# diff = ls_diff_2(2000, precision=1000)
plt.plot(diff)
plt.yscale('log')
plt.legend(['floor((3/2)^k) * (2^k -1) / (3^k - 2^k) - 1'])

# diff = ls_log_diff(1750, precision=1000)
# plt.plot(diff)
# plt.yscale('log')
# plt.legend(['log(floor((3/2)^k) * (2^k -1) / (3^k - 2^k))'])
