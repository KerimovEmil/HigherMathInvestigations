from matrix.basic_matrix import MatrixError, Matrix
from matrix.square_matrix import SquareMatrix
import copy
import itertools
import unittest


class VandermondeMatrix(Matrix):
    def __init__(self, ls_entries=None, matrix=None):
        if ls_entries is not None and isinstance(ls_entries, list):
            super().__init__(ls_entries)
        elif matrix is not None and isinstance(matrix, Matrix):
            super().__init__(ls_entries=matrix.ls_entries)
        else:
            raise MatrixError('Vandermonde matrix not instantiated correctly.')

        if not self.is_vandermonde():
            raise MatrixError('Not a vandermonde matrix.')

    def polynomial_fitting(self):
        # TODO
        pass

    def lagrange_interpolation(self):
        # TODO
        # add plotting to show the differences between vandermonde fitting and lagrange
        # https://mathrule.wordpress.com/2010/06/05/polynomial-interpolationvandermonde-matrix-lagrange-polynomial/
        pass


class SquareVandermondeMatrix(SquareMatrix, VandermondeMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

    def determinant(self):
        """
        Calculates the determinant of a Square Vandermonde matrix
        Returns:

        """
        A = copy.deepcopy(self)
        mult_res = 1
        N = self.len_row

        # the column after the column of ones
        # either the first or last column is 1s
        index_col_ones = 0 if any(x != 1 for x in A[:, -1]) else N - 1
        index_col = 1 if index_col_ones == 0 else N - 2
        ls_factors = A[:, index_col]

        for tup_comb in list(set(itertools.combinations(range(N), 2))):
            mult_res *= (ls_factors[tup_comb[1]] - ls_factors[tup_comb[0]])

        return mult_res


class TestSquareVandermondeMatrix(unittest.TestCase):
    def test_det(self):
        import numpy as np

        input = np.random.randint(8, size=6)

        A_np1 = np.vander(input, increasing=True)
        A_np2 = np.vander(input, increasing=False)

        A_test1 = SquareVandermondeMatrix(ls_entries=A_np1.tolist())
        A_test2 = SquareVandermondeMatrix(ls_entries=A_np2.tolist())

        self.assertEqual(A_test1.determinant(), np.linalg.det(A_np1))
        self.assertEqual(A_test2.determinant(), np.linalg.det(A_np2))
