import unittest


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
        return self == self.transpose()

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
        ls_new_entries = [[0 for _ in range(self.len_col)] for _ in range(self.len_row)]

        for row in range(self.len_row):
            for column in range(self.len_col):
                ls_new_entries[row][column] = self[row][column] + other[row][column]
        return self.matrix_factory(ls_new_entries)

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
        len_col = len(A.ls_entries[0])

        for i in reversed(range(1, len_col - 1)):
            reduction_fact = Matrix([[x * z[0] for x in A[i]] for z in zip(A[:i, i])])
            A[:i] = Matrix(A[:i]).__add__(reduction_fact.__neg__())
        return A

    @staticmethod
    def row_echelon_form(A):
        """
        Converts the input matrix into row echelon form.  Note, this function is recursive

        Args:
            A: <Matrix> to put into row echelon form
        Returns:
            <Matrix> row echelon form of Matrix A
        """

        len_col = len(A.ls_entries[0])

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
        if len_col > 2:
            reduction_fact = Matrix([[x * z[0][0] for x in A[0]] for z in zip(A[1:, 0:1])])
            A[1:] = Matrix(A[1:]).__add__(reduction_fact.__neg__())

        # Recursively reduce the next rows and columns
        B = Matrix.row_echelon_form(Matrix(A[1:, 1:])) if len_col > 2 else []

        # return the prior evaluated rows of A (A[:1]) appended to the result of the next set of reductions (B) plus the
        # prior columns of zero to maintain the matrix shape (A[1:, :1]
        return Matrix(A[:1] + [x + y for x, y in zip(A[1:, :1], B)])


# https://codereview.stackexchange.com/questions/233182/general-matrix-class?rq=1


class TestMatrix(unittest.TestCase):
#todo MICHELLE Add unittests for the reduced row echlon and row echlon functions
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
