import unittest
from itertools import chain
import copy


# no importing numpy


class MatrixError(Exception):
    """An exception class for Matrix"""
    pass


class Matrix:

    def __init__(self, ls_entries=None):

        # check that input is 2D
        if ls_entries and ls_entries[0]:
            if any(isinstance(x, list) for x in ls_entries[0]):  # more than 2D
                raise MatrixError('Input is not 2 dimensional')
        else:  # 1D
            raise MatrixError('Input is not 2 dimensional')

        self.ls_entries = ls_entries

        self.len_row = len(self.ls_entries)
        self.len_col = len(self.ls_entries[0])

        if not all(map(lambda x: len(x) == self.len_col, self.ls_entries)):
            raise MatrixError("Uneven columns")

        self.matrix_factory = self.get_matrix_factory()

    def get_matrix_factory(self):
        if not hasattr(self, 'matrix_factory'):
            from matrix.choose_matrix_type import MatrixFactory
            return MatrixFactory()

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return self.ls_entries[key]
        if len(key) == 2:
            row, col = key
            if isinstance(row, slice):  # for if either row only or row and col are slices
                return [x[col] for x in self.ls_entries[row]]
            elif all(isinstance(i, int) for i in key) or isinstance(col,
                                                                    slice):  # for if only col is a slice or neither row and col are slices
                return self.ls_entries[row][col]
            else:
                raise NotImplemented
        else:
            raise NotImplemented

    def __setitem__(self, index, value):
        self.ls_entries[index] = value

    def zero_matrix(self, row_dim, col_dim):
        """
        Returns a Zero matrix with row and columns
        Args:
            row_dim: <int> a positive integer representing the number of rows
            col_dim: <int> a positive integer representing the number of rows

        Returns: <Matrix> of zeros
        """
        return self.matrix_factory(ls_entries=self.zero_ls_entries(row_dim, col_dim))

    @staticmethod
    def zero_ls_entries(row_dim, col_dim):
        """
        Returns a Zero matrix with row and columns
        Args:
            row_dim: <int> a positive integer representing the number of rows
            col_dim: <int> a positive integer representing the number of rows

        Returns: <list> of zeros
        """
        assert isinstance(row_dim, int) and isinstance(col_dim, int)
        return [[0] * col_dim for _ in range(row_dim)]

    def ones_matrix(self, row_dim, col_dim):
        """
        Returns a matrix of ones with row and columns
        Args:
            row_dim: <int> a positive integer representing the number of rows
            col_dim: <int> a positive integer representing the number of rows

        Returns: <Matrix> of ones
        """
        return self.matrix_factory(ls_entries=self.ones_ls_entries(row_dim, col_dim))

    @staticmethod
    def ones_ls_entries(row_dim, col_dim):
        """
        Returns a matrix of ones with row and columns
        Args:
            row_dim: <int> a positive integer representing the number of rows
            col_dim: <int> a positive integer representing the number of rows

        Returns: <list> of ones
        """
        assert isinstance(row_dim, int) and isinstance(col_dim, int)
        return [[1] * col_dim for _ in range(row_dim)]

    def identity(self, size):  # todo see if this should be a new child class
        """
        Return an identity matrix of size n
        Args:
            size: <int> a positive integer representing the size of the identity matrix

        Returns: Identity <Matrix>
        """
        assert isinstance(size, int)
        ls_zero = self.zero_ls_entries(row_dim=size, col_dim=size)
        for i in range(size):
            ls_zero[i][i] = 1
        return self.matrix_factory(ls_entries=ls_zero)

    def transpose(self):
        return self.__class__(ls_entries=[[self[j][i] for j in range(self.len_row)] for i in range(self.len_col)])

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
            result = self.zero_ls_entries(row_dim=self.len_row, col_dim=self.len_col)

            for row in range(self.len_row):
                for col in range(self.len_col):
                    result[row][col] = self.ls_entries[row][col] * other

        elif isinstance(other, Matrix):
            if self.len_col != other.len_row:
                raise MatrixError('Multiplication dimension wrong.')
            result = self.zero_ls_entries(row_dim=self.len_row, col_dim=other.len_col)

            for i in range(self.len_row):
                # iterate through columns of Y
                for j in range(other.len_col):
                    # iterate through rows of Y
                    for k in range(other.len_row):
                        result[i][j] += self[i][k] * other[k][j]
        else:
            raise MatrixError(f'type: {type(other)} multiplication not supported.')

        return self.matrix_factory(result)

    def __eq__(self, other):
        return self.ls_entries == other.ls_entries

    def __str__(self):
        s = "\n".join([str(i) for i in [rows for rows in self.ls_entries]])
        return s

    def is_square(self):
        return self.len_row == self.len_col

    def is_symmetric(self):

        # cannot use self == self.transpose in current Matrix Factory set up or else symmetric .transpose() will overwrite
        # and all square matrices will be marked as symmetric
        return self == Matrix.transpose(Matrix(self.ls_entries))

    def is_hankel(self):
        """ Checks whether matrix is hankel.  Returns True if hankel matrix, False otherwise """
        if not self.is_symmetric():
            return False

        size = self.len_col  # since will be square matrix

        # check if the skew-diagonals are constant
        for skew_diag_idx in range(size * 2 - 1):
            skew_diag = [row[skew_diag_idx - i] for i, row in enumerate(self) if size > skew_diag_idx - i >= 0]
            first_elem = skew_diag[0]
            if not all(first_elem == elem for elem in skew_diag):
                return False

        return True

    def is_orthogonal(self):
        I = self.identity(self.len_col)
        return self.__mul__(self.transpose()) == I

    def is_vandermonde(self):
        """ Checks whether matrix is vandermonde.  Returns True if vandermonde matrix, False otherwise """
        if self.len_col == 1 or self.len_row == 1:  # uncommon use case should at least be 2x2
            return False

        ls_entries = copy.deepcopy(self.ls_entries)
        # Allow for left-to-right or right-to-left order
        if self[:, -1] == [1] * self.len_row:
            # swap column order
            ordered = list(reversed(list(zip(*ls_entries))))
            ls_entries = list([list(x) for x in zip(*ordered)])

        # Create the comparable vandermonde matrix
        vander_mat = Matrix.vander_ls_entries(input_arr=Matrix(ls_entries)[:, 1], num_cols=self.len_col)

        return ls_entries == vander_mat

    @staticmethod
    def vander_ls_entries(input_arr, num_cols=False):
        """
        Returns the ls_entries for a vandermonde matrix
        Simular functionality to np.vander - columns of the output matrix are powers of the input array

        Args:
            input_arr: <list> input array to create vandermonde matrix
            num_cols: <int> number of output columns.  Defaults to False.  If False, return square matrix
            (i.e. num_cols = len(input_arr)

        Returns:
            <list> ls_entries output for vandermonde matrix
        """

        # create a matrix of size (num_cols, len(input_arr))
        num_cols = len(input_arr) if not num_cols else num_cols

        # first column is ones
        A = Matrix([[1] * len(input_arr) for _ in range(num_cols)])
        A[1] = input_arr
        for i in range(2, num_cols):
            A[i] = [x ** i for x in input_arr]

        return Matrix(A.ls_entries).transpose().ls_entries

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

    def elem_pow(self, n):
        """
        Raise each element to the power of n

        Args:
            n: <int> to raise all elements to the power to
        Returns:
            <Matrix> raised to the power
        """
        return self.matrix_factory([list(map(lambda x: pow(x, n), row)) for row in self.ls_entries])

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self

    def __neg__(self):
        return -1 * self

    def __add__(self, other):
        assert self.len_row == other.len_row
        assert self.len_col == other.len_col
        ls_new_entries = [[0 for _ in range(self.len_col)] for _ in range(self.len_row)]

        for row in range(self.len_row):
            for column in range(self.len_col):
                ls_new_entries[row][column] = self[row][column] + other[row][column]
        return self.matrix_factory(ls_new_entries)

    def minor_matrix(self, remove_row, remove_col):
        new_matrix_array = [row[:remove_col] + row[remove_col + 1:] for row in
                            (self[:remove_row] + self[remove_row + 1:])]
        return self.matrix_factory(ls_entries=new_matrix_array)

    def diagonal(self):
        min_dim = min(self.len_row, self.len_col)
        result = []
        for i in range(min_dim):
            result.append(self.ls_entries[i][i])
        return result

    def is_zero_matrix(self):
        return all(x == 0 for x in chain(*self.ls_entries))

    @staticmethod
    def reduced_row_echelon_form(A):
        """
        Given Matrix A, return the reduced row echelon form

        Args:
            A: <Matrix> to put into row echelon form
        Returns:
            <Matrix> in reduced row echelon form
        """
        A = Matrix.row_echelon_form(A)
        len_col, len_row = len(A.ls_entries[0]), len(A.ls_entries)
        loop_from = len_col if len_row >= len_col else len_col - 1

        for i in reversed(range(1, loop_from)):
            reduction_fact = Matrix([[x * z[0] for x in A[i]] for z in zip(A[:i, i])])
            A[:i] = Matrix(A[:i]).__add__(reduction_fact.__neg__())
        return A

    @staticmethod
    def is_ref(A):
        # Check if non-zero matrix in row echelon form
        while all(x == 0 for x in A[:, 0]):
            A = Matrix(A[:, 1:])

        for i in range(len(A.ls_entries[0])):
            if A[i][i] != 1:
                return False
        return True

    @staticmethod
    def row_echelon_form(A):
        """
        Converts the input matrix into row echelon form.  Note, this function is recursive

        Args:
            A: <Matrix> to put into row echelon form
        Returns:
            <Matrix> row echelon form of Matrix A
        """
        if A.is_zero_matrix() or Matrix.is_ref(A):
            return A

        # extract index of first non-zero element in the first column
        first_col_elems = [A[i, 0] for i in range(len(A.ls_entries))]
        first_non_zero_elem_idx = next((i for i, x in enumerate(first_col_elems) if x != 0), False)
        if first_non_zero_elem_idx is False:  # there are no non-zero elements
            B = Matrix.row_echelon_form(Matrix(A[:, 1:]))  # proceed to next column

        # pivot rows if the non-zero element is in another row
        if first_non_zero_elem_idx != 0:
            pivot_row = A[first_non_zero_elem_idx].copy()
            A[first_non_zero_elem_idx], A[0] = A[0], pivot_row

        # Gaussion Elimination
        A[0] = [x / A[0][0] for x in A[0]]
        reduction_fact = [[x * z[0][0] for x in A[0]] for z in zip(A[1:, 0:1])]
        A[1:] = [list(map(lambda x, y: x - y, row_mat_1, row_mat_2)) for row_mat_1, row_mat_2 in
                 zip(A[1:], reduction_fact)]

        # Recursively reduce the next rows and columns
        if list(chain(*A[1:, 1:])):
            B = Matrix.row_echelon_form(Matrix(A[1:, 1:]))
        else:
            return A

        # return the prior evaluated rows of A (A[:1]) appended to the result of the next set of reductions (B) plus the
        # prior columns of zero to maintain the matrix shape (A[1:, :1]
        return Matrix(A[:1] + [x + y for x, y in zip(A[1:, :1], B)])

    # https://codereview.stackexchange.com/questions/233182/general-matrix-class?rq=1


class TestMatrix(unittest.TestCase):
    def test_reduced_row_echelon(self):
        import numpy as np
        from sympy import Matrix as sympy_matrix

        ls_entries_col_g_row = [[14, 0, 11, 3],
                                [22, 23, 4, 7],
                                [-12, -34, -3, -4]]
        ls_entries_square = [[16, 14, 11, 18, 11],
                             [15, 14, 8, 15, 10],
                             [15, 10, 4, 7, 15],
                             [7, 13, 9, 19, 9],
                             [12, 4, 19, 8, 9]]

        ls_entries_row_g_col = [[14, 0, 11],
                                [22, 23, 4],
                                [-12, -34, -3],
                                [4, 5, 6]]

        for ls_entries in [ls_entries_col_g_row, ls_entries_square, ls_entries_row_g_col]:
            sympy_rref, _ = sympy_matrix(ls_entries).rref()

            A = Matrix(ls_entries)
            rref = Matrix.reduced_row_echelon_form(A)

            self.assertTrue(np.allclose(rref.ls_entries, np.array(sympy_rref).astype(np.float64)))

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

        self.assertEqual(I.identity(size=3), I)

    def test_symmetric(self):
        # 3x4 matrix
        B = Matrix(ls_entries=[
            [5, 8, 1, 2],
            [6, 7, 3, 0],
            [4, 5, 9, 1]])

        self.assertFalse(B.is_symmetric())

        B = Matrix(ls_entries=[
            [5, 8, 1],
            [6, 7, 3],
            [4, 5, 9]])

        self.assertFalse(B.is_symmetric())

        B = Matrix(ls_entries=[
            [7, 6, 4],
            [6, 7, 5],
            [4, 5, 7]])

        self.assertTrue(B.is_symmetric())

    def test_hankel(self):
        B = Matrix(ls_entries=[
            [5, 8, 1],
            [6, 7, 3],
            [4, 5, 9]])

        self.assertFalse(B.is_hankel())

        B = Matrix(ls_entries=[
            [1, 2, 3, 4],
            [2, 3, 4, 0],
            [3, 4, 0, 0],
            [4, 0, 0, 0]])

        self.assertTrue(B.is_hankel())

        # Hilbert matrix
        B = Matrix(ls_entries=[
            [1, 1 / 2, 1 / 3, 1 / 4],
            [1 / 2, 1 / 3, 1 / 4, 1 / 5],
            [1 / 3, 1 / 4, 1 / 5, 1 / 6],
            [1 / 4, 1 / 5, 1 / 6, 1 / 7]])

        self.assertTrue(B.is_hankel())

    def test_orthogonal(self):
        B = Matrix([[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [1, 0, 0, 0],
                    [0, 1, 0, 0]])

        self.assertTrue(B.is_orthogonal())

        B = Matrix([[1, 0], [0, -1]])

        self.assertTrue(B.is_orthogonal())

    def test_vandermonde(self):
        B = Matrix(ls_entries=[[1, 1, 1],
                               [4, 2, 1],
                               [9, 3, 1],
                               [25, 5, 1]])

        self.assertTrue(B.is_vandermonde())

        B = Matrix(ls_entries=[[1, 1, 1, 1],
                               [8, 4, 2, 1],
                               [27, 9, 3, 1],
                               [125, 25, 5, 1]])

        self.assertTrue(B.is_vandermonde())

        B = Matrix(ls_entries=[[5, 1, 1, 1],
                               [5, 4, 2, 1],
                               [5, 9, 3, 1],
                               [5, 25, 5, 1]])

        self.assertFalse(B.is_vandermonde())
