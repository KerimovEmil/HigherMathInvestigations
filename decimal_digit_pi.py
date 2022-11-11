from special_number_series import bernoulli


def basic_factorial(x):
    """Returns the factorial of the integer x."""
    ans = 1
    while x:
        ans *= x
        x -= 1
    return ans


def helper_pi(n: int) -> float:
    m = 2*n
    num = 2 * basic_factorial(m)
    a = (1 - pow(2, -m)) * (1 - pow(3, -m)) * (1 - pow(5, -m)) * (1 - pow(7, -m))
    den = abs(float(bernoulli(m))) * a

    return pow(num / den, 1 / m) / 2


def decimal_digit_extraction(n: int) -> int:
    """Returns the nth digit of pi"""
    x = helper_pi(n-1) * pow(10, n-1)
    k = int(10 * (x - int(x)))
    return k


def f(k):
    """Bailey-Borwin-Plouffe formula"""
    k8 = 8*k
    result = 4/(k8 + 1) - 2/(k8+4) - 1/(k8+5) - 1/(k8+6)
    return result / pow(16, k)


def bbp_series(max_n=11):
    pi = 3.141592653589793
    ans = 0
    # testing convergence just for fun
    for n in range(max_n):
        ans += f(n)
    print(ans, pi, ans - pi)


if __name__ == '__main__':
    # bbp_series(max_n=11)

    pi = 3.141592653589793
    for i in range(3, 10):
        print(str(pi)[i+1], decimal_digit_extraction(i))
