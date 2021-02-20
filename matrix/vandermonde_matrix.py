from matrix.basic_matrix import MatrixError, Matrix
from matrix.square_matrix import SquareMatrix
import copy
import itertools
import unittest
from sympy import Symbol, latex, Add


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


class SquareVandermondeMatrix(SquareMatrix, VandermondeMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

    def determinant(self):
        """
        Calculates the determinant of a Square Vandermonde matrix

        Returns:
            mult_res: <float> the determinant

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

    @staticmethod
    def polynomial_fitting(pts_ls):
        """
        Create the interpolating polynomial passing through the input points using the vandermonde matrix

        Args:
            pts_ls: <list> of tuples of (x,y) points to be interpolated

        Returns:
            res: <list> of coefficients for the interpolating polynomial of degree n -1, where n is the number of input
            points, in descending order
            expr: <str> latex representation of the polynomial

        """
        # Create the vandermonde matrix from the x-coordinates of the points
        x_coords = [x[0] for x in pts_ls]
        y_coords_vector = Matrix([[x[1] for x in pts_ls]]).transpose()

        A = SquareVandermondeMatrix(Matrix.vander_ls_entries(input_arr=x_coords))
        B = A.matrix_factory([list(x) for x in [*map(lambda rows: itertools.chain(*rows), zip(*[A, y_coords_vector]))]])

        # solve the linear system now, Ax = b
        res = Matrix.reduced_row_echelon_form(B)[:, -1]
        res.reverse()

        poly_exp_ls = list(reversed(range(len(x_coords))))
        x = Symbol("x")

        # create expression
        expr = Add(*[a * x ** b for a, b in zip(res, poly_exp_ls)])

        return res, latex(expr)


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

    def test_interpolating_poly(self):
        input_ls = [(1, -6), (2, 2), (4, 12)]
        res_test, expr_test = SquareVandermondeMatrix.polynomial_fitting(input_ls)

        res_expected = [-1.0, 11.0, -16.0]
        expr_expected = '- 1.0 x^{2} + 11.0 x - 16.0'
        self.assertEqual(res_test, res_expected)
        self.assertEqual(expr_test, expr_expected)
