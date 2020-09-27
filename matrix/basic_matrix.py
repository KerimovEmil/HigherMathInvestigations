import unittest


class Matrix:
    def __init__(self, ls_entries=None):
        self.ls_entries = ls_entries
        # todo add tuple definition
        # todo add validation of dimension

        self.len_row = len(self.ls_entries)
        self.len_col = len(self.ls_entries[0])

    def transpose(self):
        return Matrix(ls_entries=[[self[j][i] for j in range(self.len_row)] for i in range(self.len_col)])

    def __getitem__(self, key):
        return self.ls_entries[key]

    def __mul__(self, other):
        # self * other
        result = [[0 for _ in range(other.len_col)] for _ in range(self.len_row)]
        for i in range(self.len_row):
            # iterate through columns of Y
            for j in range(other.len_col):
                # iterate through rows of Y
                for k in range(other.len_row):
                    result[i][j] += self[i][k] * other[k][j]
        return Matrix(result)

    def __eq__(self, other):
        return self.ls_entries == other.ls_entries

    def __str__(self):
        s = "\n".join([str(i) for i in [rows for rows in self.ls_entries]])
        return s

    def is_square(self):
        return self.len_row == self.len_col

    def __mod__(self, mod):
        if mod:
            for i in range(len(self.ls_entries)):
                for j in range(len(self.ls_entries[0])):
                    self.ls_entries[i][j] %= mod
        return self

    def __pow__(self, n, mod=None):
        assert (n > 0)
        if n == 1:
            return self.__mod__(mod)
        half = self.__pow__(n >> 1, mod)
        if n & 1 == 1:  # if odd
            return half.__mul__(half).__mul__(self).__mod__(mod)
        else:  # if even
            return half.__mul__(half).__mod__(mod)

    def det(self):
        raise NotImplementedError

    def square_root(self):
        raise NotImplementedError

    def LU_decomposition(self):
        raise NotImplementedError

    def eigenvectors(self):
        raise NotImplementedError

    def __invert__(self):
        raise NotImplementedError

    def diagonal(self):
        raise NotImplementedError

# https://codereview.stackexchange.com/questions/233182/general-matrix-class?rq=1


class TestMatrix(unittest.TestCase):
    def test_simple_multiplication(self):
        A = Matrix(ls_entries=[[1, 2], [1, 3]])
        B = Matrix(ls_entries=[[1, 0], [0, 1]])

        self.assertEqual(A*B, A)

    def test_non_equal_multiplication(self):
        # 3x3 matrix
        A = Matrix(ls_entries=[
            [12, 7, 3],
            [4, 5, 6],
            [7, 8, 9]])
        # 3x4 matrix
        B = Matrix(ls_entries=[
            [5, 8, 1, 2],
            [6, 7, 3, 0],
            [4, 5, 9, 1]])
        # 3x4 matrix
        C = Matrix(
            ls_entries=[
                [114, 160, 60, 27],
                [74, 97, 73, 14],
                [119, 157, 112, 23]])

        self.assertEqual(A*B, C)

    def test_transpose(self):
        # 3x4 matrix
        B = Matrix(ls_entries=[
            [5, 8, 1, 2],
            [6, 7, 3, 0],
            [4, 5, 9, 1]])

        B_t = Matrix(ls_entries=
                     [[5, 6, 4],
                      [8, 7, 5],
                      [1, 3, 9],
                      [2, 0, 1]])
        self.assertEqual(B.transpose(), B_t)
        # test double transpose
        self.assertEqual(B.transpose().transpose(), B)
