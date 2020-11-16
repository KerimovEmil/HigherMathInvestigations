from matrix.basic_matrix import Matrix, MatrixError
import unittest
import copy
from numpy import roots as numpy_roots
import itertools


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
            principal = self.matrix_factory(
                [[A[j][x] for x in range(len(A[j])) if x != 0] for j in range(len(A.ls_entries)) if
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
        """

        DELETE AFTER:
        np.linalg.eig(np_arr)
(array([ 58.83444011,   7.90991199,   3.85703574,   1.61848207,
       -10.2198699 ]), array([[-0.5353104 ,  0.30931754,  0.30174974, -0.04065062, -0.03246036],
       [-0.47786848,  0.41176689,  0.66749564,  0.78345659,  0.09460916],
       [-0.39347121, -0.32543656, -0.02559092,  0.18621625,  0.71601608],
       [-0.4260566 ,  0.207349  , -0.3509878 , -0.52917708, -0.0424597 ],
       [-0.38566912, -0.76542022, -0.58270862, -0.26426962, -0.68957492]]))

58.83444010710231 [1.0, 0.0, 0.0, 0.0, -1.3880043094064085, -0.0]
[0.0, 1.0, 0.0, 0.0, -1.239063382838598, -0.0]
[0.0, 0.0, 1.0, 0.0, -1.0202300265048931, -0.0]
[0.0, 0.0, 0.0, 1.0, -1.104720544909695, -0.0]
[0.0, 0.0, 0.0, 0.0, 1.0, -0.0]
        :return:
        """
        eig_vals_ls = self.eigenvalues()
        A = self.matrix_factory(copy.deepcopy(self.ls_entries))
        I = A.identity(size=A.size)
        eigs_result = []

        for eig_val in eig_vals_ls:
            B = A.__add__(eig_val * I.__neg__())
            # Equate to 0
            zero_vector = self.zero_matrix(self.size, 1)
            # add zero vector as end column
            B = self.matrix_factory(
                [list(x) for x in [*map(lambda rows: itertools.chain(*rows), zip(*[B, zero_vector]))]])
            ref = Matrix.row_echelon_form(B)

            # partially reduce further, to second last col
            for i in reversed(range(1, A.size - 1)):
                reduction_fact = Matrix([[x * z[0] for x in ref[i]] for z in zip(ref[:i, i])])
                ref[:i] = Matrix(ref[:i]).__add__(reduction_fact.__neg__())

            eig_vec = ref[:, -2]
            # normalize
            norm_fact = (sum([x ** 2 for x in eig_vec])) ** 0.5
            eig_vec = [x / norm_fact for x in eig_vec]

            # last element is opposite sign
            eig_vec[-1] *= -1
            eigs_result.append(eig_vec)
        return eig_vals_ls, self.matrix_factory(ls_entries=eigs_result).transpose()

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

    def test_eigenvalues_eigenvectors(self):
        ls_entries = [[16, 14, 11, 18, 11],
                      [15, 14, 8, 15, 10],
                      [15, 10, 4, 7, 15],
                      [7, 13, 9, 19, 9],
                      [12, 4, 19, 8, 9]]

        from numpy.linalg import eig as numpy_eig
        from numpy import array as numpy_array, allclose

        A = SquareMatrix(ls_entries)
        A_numpy = numpy_array(ls_entries)

        eigs, eig_vects = A.eigenvectors()
        np_eigs, np_eig_vects = numpy_eig(A_numpy)

        eig_vects = sorted(eig_vects.transpose().ls_entries, key=lambda x: x[0])
        np_eig_vects = sorted(np_eig_vects.T, key=lambda x: x[0])

        # Just need to test for the eigenvectors, if those match then the eigenvalues would've also been the same
        # Tests to the default tolerances of rtol=1e-05, atol=1e-08
        self.assertTrue(allclose(eig_vects, np_eig_vects))
