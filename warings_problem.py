# if 2^k frac( (3/2)^k ) + floor( (3/2)^k ) <= 2^k for all k then
# g(k) = 2^k + floor( (3/2)^k ) - 2

# this eq has been proved valid for k < 471,600,000

# 2^k ((3/2)^k - floor( (3/2)^k )) + floor( (3/2)^k ) <= 2^k
# 3^k + (1-2^k) * floor( (3/2)^k )) <= 2^k
# 3^k - 2^k <= (2^k - 1) * floor( (3/2)^k ))

# (3/2)^k - 1 <= (1 - 2^-k) * floor( (3/2)^k ))

# floor( (3/2)^k )) => ((3/2)^k - 1) / (1 - 2^-k)
# floor( (3/2)^k )) / (3/2)^k => (1- (2/3)^k) / (1 - 2^-k)
# floor((3/2)^k )) / (3/2)^k => (6^k- 4^k) / (6^k - 3^k)


import matplotlib.pyplot as plt
import math

m = 100

# plotting the difference
# diff = [int((3/2)**i) * (2/3)**i - (6**i - 4**i) / (6**i - 3**i) for i in range(2, m)]
diff = [int(math.pow(1.5, i)) * math.pow(2/3, i)
        -
        (1 - math.pow(2/3, i)) / (1 - math.pow(2, -i))
        for i in range(2, m)]
plt.plot(diff)
plt.yscale('log')
plt.ylabel('Inequality greater than 0')
plt.xlabel('k')
plt.legend(['floor((3/2)^k )) / (3/2)^k - (6^k- 4^k) / (6^k - 3^k)'])

# plotting the floor((3/2)^k )) / (3/2)^k
floor_ratio = [1 - int((3/2)**i) / (3/2)**i for i in range(2, m)]
plt.plot(floor_ratio)
plt.yscale('log')
plt.ylabel('Behaviour of floor((3/2)^k )) / (3/2)^k')
plt.xlabel('k')
plt.legend(['1 - floor((3/2)^k )) / (3/2)^k '])

# plotting the (6^k- 4^k) / (6^k - 3^k)
power_ratio = [1 - (6**i - 4**i) / (6**i - 3**i) for i in range(2, m)]
plt.plot(power_ratio)
plt.yscale('log')
plt.ylabel('Behaviour of (6^k- 4^k) / (6^k - 3^k)')
plt.xlabel('k')
plt.legend(['1 - (6^k- 4^k) / (6^k - 3^k)'])

from decimal import *
a = 3
b = a - 1
m = 1750


def diffa(a, b):
    getcontext().prec = 1000
    a = Decimal(a)
    b = Decimal(b)
    return [int(a**i / b**i) * (b**i - 1) / (a**i - b**i) - 1 for i in range(2, m)]


diff = diffa(a, b)
plt.plot(diff)
plt.yscale('log')
