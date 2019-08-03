import fractions
# B_m = sum_{k=0 to m} (m choose k) B_k


def combin(n, r):
    """A fast way to calculate binomial coefficients by Andrew Dalke (contrib)."""
    if 0 <= r <= n:
        ntok = 1
        rtok = 1
        for t in range(1, min(r, n - r) + 1):
            ntok *= n
            rtok *= t
            n -= 1
        return ntok // rtok  # bit-wise operation
    else:
        return 0


class BernoulliNumber:
    def __init__(self):
        self.cache_values = {0: fractions.Fraction(1, 1)}  # B_0 = 1

    def get(self, n):
        if n in self.cache_values:
            return self.cache_values[n]

        if n % 2 != 0 and n != 1:
            self.cache_values[n] = fractions.Fraction(0, 1)
            return fractions.Fraction(0, 1)

        b_sum = 0
        for k in range(n):
            # b_sum += combin(n + 1, k) * self.get(k)
            choose_k = fractions.Fraction(combin(n+1, k))
            b_sum += choose_k * self.get(k)
        # b_sum *= -1/(n+1)
        b_sum *= fractions.Fraction(-1, n+1)
        self.cache_values[n] = b_sum
        return b_sum


if __name__ == '__main__':
    b = BernoulliNumber()
    assert b.get(0) == 1
    assert b.get(20) == fractions.Fraction(-174611, 330)
    assert b.get(40) == fractions.Fraction(-261082718496449122051, 13530)

    for i in range(21):
        print(i, b.get(i))
