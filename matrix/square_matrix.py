from matrix.basic_matrix import Matrix, MatrixError
import unittest
import copy


class SquareMatrix(Matrix):
    def __init__(self, ls_entries=None, matrix=None):
        """
        Construct a Square Matrix class
        Args:
            ls_entries: given a list of lists construct a matrix
            matrix: given a Matrix type
        """
        if ls_entries is not None and isinstance(ls_entries, list):
            super().__init__(ls_entries)
        elif matrix is not None and isinstance(matrix, Matrix):
            super().__init__(ls_entries=matrix.ls_entries)
        else:
            raise MatrixError('Square matrix not instantiated correctly.')

        if not self.is_square():
            raise MatrixError('Not a square matrix.')
        self.size = self.len_row

    def determinant(self):
        """Computes the determinate of the square matrix"""
        # base case for 2x2 matrix
        if self.size == 2:
            # a*d - b*c
            return self[0][0] * self[1][1] - self[0][1] * self[1][0]

        determinant = 0
        for col in range(self.size):
            determinant += ((-1) ** col) * self[0][col] * self.minor_matrix(0, col).determinant()

        return determinant

    def inverse(self):
        determinant = self.determinant()
        if determinant == 0:
            raise MatrixError("Matrix is singular")

        # special case for 2x2 matrix:
        if self.size == 2:
            return [[self[1][1] / determinant, -1 * self[0][1] / determinant],
                    [-1 * self[1][0] / determinant, self[0][0] / determinant]]

        # find matrix of minors
        minors = self.zero_matrix(row_dim=self.size, col_dim=self.size)

        for row in range(self.size):
            for col in range(self.size):
                minor = self.minor_matrix(row, col)
                minors[row][col] = ((-1) ** (row + col)) * minor.determinant()

        t_minors = minors.transpose()

        for row in range(t_minors.len_row):
            for col in range(t_minors.len_row):
                t_minors[row][col] /= determinant

        return t_minors

    def square_root(self):
        raise NotImplementedError

    def LU_decomposition(self):
        """
        Returns: the LU decomposition of the self.ls_entries matrix
        self.ls_entries matrix must be square
            L, U where [A] = [L][U]
            L: <Matrix> m x m lower triangular matrix
            U: <Matrix> m x m upper triangular matrix
        """

        # L and U are both square matrices
        L = self.identity(self.len_row).ls_entries
        U = [[0] * self.len_row for x in range(self.len_row)]
        A = SquareMatrix(copy.deepcopy(self.ls_entries))

        # Need to pivot first for stability
        for i in range(self.len_row):
            minor_col_i = list(map(abs, A[i:self.len_row, i]))
            max_row_index = minor_col_i.index(max(minor_col_i)) + i
            A[i], A[max_row_index] = A[max_row_index], A[i]

        for i in range(self.len_row):
            for j in range(i, self.len_row):  # update matrix after pivoting
                sum_upper = sum([L[i][k] * U[k][j] for k in range(i)])
                sum_lower = sum([L[j][k] * U[k][i] for k in range(i)])
                U[i][j] = A[i][j] - sum_upper
                L[j][i] = (A[j][i] - sum_lower) / U[i][i] if U[i][i] != 0 else 0

        return SquareMatrix(L), SquareMatrix(U)

    def get_toeplitz_matrix_berkowitz(self, a_0_0, row_vector, col_vector, principal):
        """

        :return: (n + 1) x n Toeplitz lower triangular matrix
        """
        size = len(principal.ls_entries) + 1  # size of the matrix that the row, vector, and principal were extracted
        L = Matrix(self.identity(self.size).ls_entries + [[0] * size])

        principal_pows = [principal.__pow__(x) for x in range(1, size - 1)]
        print(len(principal.ls_entries))
        # first and second diagonals will just be -a11 and the -RS (no principal)
        # then the ones after will use the principal pows times the -RS
        # then do the loop to put values in, might be able to do the i-j is a certain value concept, not that sure yet



        print(principal_pows)
        diag_vals = []
        for i in range(1, self.size + 1):  # rows
            for j in range(i):
                if i - j == 1:
                    L[i][j] = -a_0_0
                else:
                    print('test')
                    print(principal.__pow__(2))
                    prince_test = principal.__pow__(2)
                    print(row_vector.__mul__(prince_test))
                    # L[i][j] = -row_vector.__mul__(principal).__pow__(diag_counter - 2).__mul__(col_vector)

        return L

    def char_eqn_berkowitz(self):
        # https://en.wikipedia.org/wiki/Samuelson%E2%80%93Berkowitz_algorithm
        # todo
        # probably try to get the characteristic eqm using this algorithm
        # and then use some root finding method to get the eigenvalues...
        # not sure yet though
        A = SquareMatrix(copy.deepcopy(self.ls_entries))
        for i in range(self.size):
            # row_vector = Matrix([A[i][x] for x in range(1, self.size)])
            row_vector = Matrix([[2, 1]])
            col_vector = Matrix([A[:, i]])
            principal = Matrix([[A[j][x] for x in range(len(A[j])) if x != i] for j in range(self.size) if
                                j != i])  # todo make this cleaner
            a_elem = A[0][0]
            C = self.get_toeplitz_matrix_berkowitz(a_elem, row_vector, col_vector, principal)
        return 0

    def eigenvalues(self):
        """
        Uses power iteration to calculate the eigenvalues

        :return:
        """

        # get the characteristic equation

        # power iteration?
        # inverse iteration?
        # figure out how to solve characteristic equation without using numpy

        pass

    def eigenvectors(self):
        raise NotImplementedError


class TestSquareMatrix(unittest.TestCase):
    def test_det(self):
        A = SquareMatrix(ls_entries=[[1, 2, 1], [1, 3, 2], [1, 5, 2]])

        self.assertEqual(A.determinant(), -2)

    def test_inverse(self):
        A = SquareMatrix(ls_entries=[[1, 2, 1], [1, 3, 2], [1, 5, 2]])

        A_I = 0.5 * SquareMatrix(ls_entries=[[4, -1, -1], [0, -1, 1], [-2, 3, -1]])
        self.assertEqual(A.inverse(), A_I)

    def test_lu_decomposition(self):
        A = SquareMatrix(ls_entries=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        L, U = A.LU_decomposition()
        multiply = L.__mul__(U)
        self.assertEqual(set(map(tuple, A.ls_entries)), set(map(tuple, multiply)))
