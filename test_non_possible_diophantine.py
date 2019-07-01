# Quickly see if some diophantine equations have no solution by looking at certain simple modules

# Example equation: 15 x^2 - 7 y^2 = 1
# A counter example is found using mod 5
# Proof:
# 15 mod 5 = 0, -7 mod 5 = 3
# 3 y^2 == 1 mod 5
# Looking at all squares mod 5 we get:
# 3 * 0^2 = 3 * 0 == 0 mod 5 != 1
# 3 * 1^2 = 3 * 1 == 3 mod 5 != 1
# 3 * 2^2 = 3 * 4 == 2 mod 5 != 1
# 3 * 3^2 = 3 * 9 == 2 mod 5 != 1
# 3 * 4^2 = 3 * 16 == 3 mod 5 != 1
# therefore no solutions can exist by looking at mod 5


class IntVariable:
    def __init__(self, mod, power=None, coeff=None, possible_values=None):
        self.mod = mod
        if power is not None:
            self.power = power
            self.coeff = coeff
            self.possible_values = self.calc_vals()
        if possible_values is not None:
            self.possible_values = possible_values

    def calc_vals(self):
        return set([(self.coeff * (i ** self.power)) % self.mod for i in range(self.mod)])

    def values(self):
        return self.possible_values

    def __mul__(self, other):
        # todo implement, remove coeff
        pass

    def __add__(self, other):
        if self.mod != other.mod:
            raise NotImplementedError("Addition not defined for different modules.")
        ans = set()
        for x in self.possible_values:
            for y in other.possible_values:
                ans.add((x+y) % self.mod)
        return IntVariable(mod=self.mod, possible_values=ans)

    def __repr__(self):
        return str(self.possible_values)


class IntPolynomial:
    def __init__(self, ls_coeff, ls_power):
        assert len(ls_coeff) == len(ls_power)
        self.ls_coeff = ls_coeff
        self.ls_power = ls_power

        self.poly_possible = None

    def test(self, mod, constant_value=0):
        ls_v = [IntVariable(mod=mod, power=p, coeff=c) for c, p in zip(self.ls_coeff, self.ls_power)]

        # todo replace with a prettier sum
        self.poly_possible = ls_v[0]
        for i in range(len(ls_v) - 1):
            self.poly_possible += ls_v[i + 1]

        return constant_value in self.poly_possible.values()


if __name__ == "__main__":
    # 15 x^2 - 7 y^2 - 1 = 0
    ls_coeff = [15, -7, -1]
    ls_powers = [2, 2, 0]
    for m in range(2, 10):
        mod_m_possible = IntPolynomial(ls_coeff, ls_powers).test(m)
        if mod_m_possible is False:
            print("Polynomial fails at mod: {}".format(m))
            break

    # 15 x^2 - 7 y^2 = 1
    ls_coeff = [15, -7]
    ls_powers = [2, 2]
    for m in range(2, 10):
        mod_m_possible = IntPolynomial(ls_coeff, ls_powers).test(m, 1)
        if mod_m_possible is False:
            print("Polynomial fails at mod: {}".format(m))
            break

    # x^2 − 3y^2 - 175 = 0
    ls_coeff = [1, -3, -175]
    ls_powers = [2, 2, 0]
    for m in range(2, 10):
        mod_m_possible = IntPolynomial(ls_coeff, ls_powers).test(m)
        if mod_m_possible is False:
            print("Polynomial fails at mod: {}".format(m))
            break