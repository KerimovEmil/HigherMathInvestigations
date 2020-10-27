from matrix.basic_matrix import Matrix, MatrixError
import unittest
import copy
from numpy import roots as numpy_roots


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
        ls_minors = self.zero_ls_entries(row_dim=self.size, col_dim=self.size)

        for row in range(self.size):
            for col in range(self.size):
                minor = self.minor_matrix(row, col)
                ls_minors[row][col] = ((-1) ** (row + col)) * minor.determinant()

        t_minors = self.matrix_factory(ls_entries=ls_minors).transpose()

        for row in range(t_minors.len_row):
            for col in range(t_minors.len_row):
                t_minors[row][col] /= determinant

        return t_minors

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
        A = self.matrix_factory(copy.deepcopy(self.ls_entries))

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

        return self.matrix_factory(L), self.matrix_factory(U)

    def get_toeplitz_matrix_berkowitz(self, a_0_0, row_vector, col_vector, principal, col_num):
        """
        Gets the (n+1) x n Toeplitz matrix associated with a n x n matrix, already partitioned into its first element,
        row and column vectors, and principal submatrix

        The n x n vector is partitioned per below:
        A = |a_0_0     |row_vector|
            |col_vector|principal |

        Args:
            a_0_0: <int> or <float> element A[0][0] of the n x n matrix
            row_vector: <Matrix> row vector of the n x n  matrix
            col_vector: <Matrix> column vector of the n x n matrix
            principal: <Matrix> principal submatrix of the n x n matrix
            col_num: <int> number of columns in the Toeplitz matrix, equal to n
        Returns:
            L: <Matrix> (n + 1) x n Toeplitz lower triangular matrix

        For more info, see:
        https://handwiki.org/wiki/Samuelson%E2%80%93Berkowitz_algorithm#:~:text=In%20mathematics%2C%20the%20Samuelson%E2%80%93Berkowitz,commutative%20ring%20without%20zero%20divisors.
        """
        row_num = col_num + 1  # number of rows of Toeplitz matrix
        L = self.matrix_factory(self.identity(col_num).ls_entries + [[0] * (col_num)])
        products = [-a_0_0, -row_vector.__mul__(col_vector.transpose())[0][0]]

        # get the diagonal entries where taking the power of the principal submatrix is required
        if row_num > 3:
            principal_pows = [principal] + [principal.__pow__(x) for x in range(2, col_num)]
            products += [-row_vector.__mul__(x).__mul__(col_vector.transpose())[0][0] for x in principal_pows]

        # Assign values to the diagonals
        for i in range(1, row_num):
            for j in range(col_num):
                if i - j > 0:
                    L[i][j] = products[(i - j) - 1]
        return L

    def char_eqn_berkowitz(self):
        """
        Gets the coefficients of the characteristic polynomial using the Berkowitz Algorithm
        # https://en.wikipedia.org/wiki/Samuelson%E2%80%93Berkowitz_algorithm


        Returns:
            <Matrix> vector containing the coefficients of the characteristic polynomial from highest to lowest degree
        """
        A = self.matrix_factory(copy.deepcopy(self.ls_entries))
        C_ls = []

        for i in range(self.size - 1):
            # partition the n x n matrix accordingly
            row_vector = self.matrix_factory([A[0, 1:]])
            col_vector = self.matrix_factory([A[1:, 0]])
            principal = self.matrix_factory([[A[j][x] for x in range(len(A[j])) if x != 0] for j in range(len(A.ls_entries)) if
                                j != 0])
            a_elem = A[0][0]
            C = self.get_toeplitz_matrix_berkowitz(a_elem, row_vector, col_vector, principal, col_num=len(A.ls_entries))
            C_ls.append(C)
            A = principal  # next Toeplitz matrix is found for subsequent principal submatrix

        # last one is just [1,-a,1,1]
        C_ls.append(self.matrix_factory(ls_entries=[[1, -A[0][0]]]).transpose())

        # create the resulting vector that stores the coefficients of the characteristic polynomial
        result = C_ls[0]
        C_ls.pop(0)
        while C_ls:
            result = result.__mul__(C_ls[0])
            C_ls.pop(0)

        return result.transpose()

    def eigenvalues(self):
        """
        Uses the numpy roots function to calculate the roots of the characteristic polynomial
        Characteristic equation is found using the Berkowitz Algorithm
        """
        char_eqn = self.char_eqn_berkowitz()
        return list(numpy_roots(char_eqn.ls_entries[0]))  # this uses the numpy roots function

    def eigenvectors(self):
        raise NotImplementedError

    def diagonalize(self):
        raise NotImplementedError

    def square_root(self):
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

    def test_eigenvalues(self):
        ls_entries = [[16, 14, 11, 18, 11],
                      [15, 14, 8, 15, 10],
                      [15, 10, 4, 7, 15],
                      [7, 13, 9, 19, 9],
                      [12, 4, 19, 8, 9]]

        from numpy.linalg import eig as numpy_eig
        from numpy import array as numpy_array

        A = SquareMatrix(ls_entries)
        A_numpy = numpy_array(ls_entries)

        eigs = A.eigenvalues()
        np_eigs, np_eig_vects = numpy_eig(A_numpy)

        self.assertTrue(all([abs(i - j) <= 1e-5 for i, j in zip(sorted(list(np_eigs)), sorted(eigs))]))
