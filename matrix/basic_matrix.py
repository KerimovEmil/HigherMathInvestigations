import unittest
import numpy as np  # do not use any numpy matrix functions
import copy


class MatrixError(Exception):
    """An exception class for Matrix"""
    pass


class Matrix:

    def __init__(self, ls_entries=None):

        # check that input is 2D
        if len(np.array(ls_entries).shape) != 2:
            raise MatrixError('Input is not 2 dimensional')

        self.ls_entries = ls_entries

        self.len_row = len(self.ls_entries)
        self.len_col = len(self.ls_entries[0])

        if not all(map(lambda x: len(x) == self.len_col, self.ls_entries)):
            raise MatrixError("Uneven columns")

    @staticmethod
    def zero_matrix(row_dim, col_dim):
        """
        Returns a Zero matrix with row and columns
        Args:
            row_dim: <int> a positive integer representing the number of rows
            col_dim: <int> a positive integer representing the number of rows

        Returns: <Matrix> of zeros
        """
        assert isinstance(row_dim, int) and isinstance(col_dim, int)
        ls_entries = [[0] * col_dim for _ in range(row_dim)]
        return Matrix(ls_entries=ls_entries)

    @staticmethod
    def identity(size):
        """
        Return an identity matrix of size n
        Args:
            size: <int> a positive integer representing the size of the identity matrix

        Returns: Identity <Matrix>
        """
        assert isinstance(size, int)
        zero = Matrix.zero_matrix(row_dim=size, col_dim=size)
        for i in range(size):
            zero[i][i] = 1
        return __class__(ls_entries=zero.ls_entries)

    def transpose(self):
        return self.__class__(ls_entries=[[self[j][i] for j in range(self.len_row)] for i in range(self.len_col)])

    def __getitem__(self, key):
        return self.ls_entries[key]

    def __rmul__(self, other):
        # other * self
        if isinstance(other, (int, float, complex)):
            return self * other
        elif isinstance(other, Matrix):
            # A * B = (B^T * A^T)^T
            return (self.transpose() * other.transpose()).transpose()

    def __mul__(self, other):
        # self * other

        if isinstance(other, (int, float, complex)):
            result = [[0 for _ in range(self.len_col)] for _ in range(self.len_row)]

            for row in range(self.len_row):
                for col in range(self.len_col):
                    result[row][col] = self.ls_entries[row][col] * other

        elif isinstance(other, Matrix):
            if self.len_col != other.len_row:
                raise MatrixError('Multiplication dimension wrong.')

            result = [[0 for _ in range(other.len_col)] for _ in range(self.len_row)]
            for i in range(self.len_row):
                # iterate through columns of Y
                for j in range(other.len_col):
                    # iterate through rows of Y
                    for k in range(other.len_row):
                        result[i][j] += self[i][k] * other[k][j]
        else:
            raise MatrixError(f'type: {type(other)} multiplication not supported.')

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

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self

    def __neg__(self):
        return -1 * self

    def __add__(self, other):
        assert self.len_row == other.len_row
        assert self.len_col == other.len_col
        new_matrix = self.zero_matrix(self.len_row, self.len_col)

        for row in range(self.len_row):
            for column in range(self.len_col):
                new_matrix[row][column] = self[row][column] + other[row][column]
        return new_matrix

    def minor_matrix(self, remove_row, remove_col):
        new_matrix_array = [row[:remove_col] + row[remove_col + 1:] for row in
                            (self[:remove_row] + self[remove_row + 1:])]
        return self.__class__(ls_entries=new_matrix_array)

    def diagonal(self):
        min_dim = min(self.len_row, self.len_col)
        result = []
        for i in range(min_dim):
            result.append(self.ls_entries[i][i])
        return result

    def LU_decomp(self):
        """
        Returns: the LU decomposition of the self.ls_entries matrix
            L, U where [A] = [L][U]
            L: <Matrix> m x m lower triangular matrix
            U: <Matrix> upper triangular matrix
        """
        L = self.identity(self.len_row).ls_entries  # L is square matrix
        U = [[0] * self.len_col for x in range(self.len_row)]
        A = copy.deepcopy(self.ls_entries)

        # Need to pivot first for stability
        for i in range(self.len_row):
            max_row_index = np.argmax(abs(np.asarray(A)[i:self.len_row, i])) + i #todo check with Emil if use of np.asarray is allowed
            A[i], A[max_row_index] = A[max_row_index], A[i]

        for i in range(self.len_row):
            for j in range(i, self.len_row):  # update matrix after pivoting
                sum_upper = sum([L[i][k] * U[k][j] for k in range(i)])
                sum_lower = sum([L[j][k] * U[k][i] for k in range(i)])
                U[i][j] = A[i][j] - sum_upper
                L[j][i] = (A[j][i] - sum_lower) / U[i][i] if U[i][i] != 0 else 0

        return Matrix(L), Matrix(U)


# https://codereview.stackexchange.com/questions/233182/general-matrix-class?rq=1


class TestMatrix(unittest.TestCase):
    def test_lu_decomposition(self):
        A = Matrix(ls_entries=[
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]])

        L, U = A.LU_decomp()
        multiply = L.__mul__(U)
        self.assertEqual(set(map(tuple, A.ls_entries)), set(map(tuple, multiply)))

    def test_simple_multiplication(self):
        A = Matrix(ls_entries=[[1, 2], [1, 3]])
        B = Matrix(ls_entries=[[1, 0], [0, 1]])

        self.assertEqual(A * B, A)

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

        self.assertEqual(A * B, C)

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

    def test_add(self):
        # 3x4 matrix
        B = Matrix(ls_entries=[
            [5, 8, 1, 2],
            [6, 7, 3, 0],
            [4, 5, 9, 1]])

        self.assertEqual(B + B, 2 * B)

    def test_identity(self):
        I = Matrix(ls_entries=[
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

        self.assertEqual(Matrix.identity(size=3), I)
