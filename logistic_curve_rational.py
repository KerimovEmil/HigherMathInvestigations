import fractions

# x_{n+1} = lam * x_n * (1-x_n)


if __name__ == '__main__':

    lam = fractions.Fraction(30, 10)
    f = lambda x: lam * x * (1 - x)

    # iterations = 100
    iterations = 10
    x_i = fractions.Fraction(1, 2)

    for i in range(iterations):
        print(x_i)
        x_i = f(x_i)
    print(x_i)

