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
    def __init__(self, mod, possible_values=None):
        self.mod = mod
        self.possible_values = possible_values if possible_values is not None else set(range(self.mod))

    def __mul__(self, other):
        if isinstance(other, int):
            other_list = [other]
        elif isinstance(other, IntVariable):
            if self.mod != other.mod:
                raise NotImplementedError("Addition not defined for different modules.")
            other_list = other.possible_values
        else:
            raise NotImplementedError("For type: {}".format(type(other)))

        possible_values = {(x * y) % self.mod for x in self.possible_values for y in other_list}
        return IntVariable(mod=self.mod, possible_values=possible_values)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, power):
        possible_values = {(i ** power) % self.mod for i in self.possible_values}
        return IntVariable(mod=self.mod, possible_values=possible_values)

    def __add__(self, other):
        if isinstance(other, int):
            other_list = [other]
        elif isinstance(other, IntVariable):
            if self.mod != other.mod:
                raise NotImplementedError("Addition not defined for different modules.")
            other_list = other.possible_values
        else:
            raise NotImplementedError("For type: {}".format(type(other)))

        possible_values = {(x + y) % self.mod for x in self.possible_values for y in other_list}
        return IntVariable(mod=self.mod, possible_values=possible_values)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-1*other)

    def __rsub__(self, other):
        return -1 * self.__sub__(other)

    def __repr__(self):
        return str(self.possible_values)


class IntPolynomial:
    def __init__(self, ls_coeff, ls_power):
        assert len(ls_coeff) == len(ls_power)
        self.ls_coeff = ls_coeff
        self.ls_power = ls_power

        self.poly_possible = None

    def test(self, mod, constant_value=0):
        ls_v = [c*IntVariable(mod=mod)**p for c, p in zip(self.ls_coeff, self.ls_power)]
        self.poly_possible = sum(ls_v)
        return constant_value in self.poly_possible.possible_values


def break_mod_loop(ls_coeff, ls_powers, max_mod=10):
    """Only works for independent terms"""
    for m in range(2, max_mod):
        mod_m_possible = IntPolynomial(ls_coeff, ls_powers).test(m)
        if mod_m_possible is False:
            # print("Polynomial fails at mod: {}".format(m))
            return m
    return None


if __name__ == "__main__":
    # 15 x^2 - 7 y^2 - 1 = 0
    assert break_mod_loop(ls_coeff=[15, -7, -1], ls_powers=[2, 2, 0]) == 3

    # x^2 − 3y^2 - 175 = 0
    assert break_mod_loop(ls_coeff=[1, -3, -175], ls_powers=[2, 2, 0]) == 4

    # x^2 - 15 y^2 - 22 = 0
    assert break_mod_loop(ls_coeff=[1, -15, -22], ls_powers=[2, 2, 0]) == 5

    # x^2 + 4x + 1 - 4y^2 = 0
    assert break_mod_loop(ls_coeff=[1, 4, 1, -4], ls_powers=[2, 1, 0, 2]) == 4

    # easiest implementation and most accurate
    break_mod = None
    for m in range(2, 1000):
        poly_vals = {(x ** 3 - 3 * x * (y ** 2) + y ** 3 - 2891) % m for x in range(m) for y in range(m)}
        if 0 not in poly_vals:
            break_mod = m
            break
    assert break_mod == 9
