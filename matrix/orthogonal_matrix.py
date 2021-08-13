from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError, Matrix
import unittest


class OrthogonalMatrix(SquareMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

        if not self.is_orthogonal():
            raise MatrixError('Not an orthogonal matrix.')

    def inverse(self):
        """ Inverse is the same as the transpose """
        return self.transpose()

    def determinant(self):
        return [1, -1]

    def solve_linear_system(self, b):
        """
        Args:
            b: <matrix.basic_matrix.Matrix> constant matrix

        Returns:
            x: <matrix.basic_matrix.Matrix> vector solution of linear system
        """
        return self.transpose().__mul__(b)


class TestOrthogonalMatrix(unittest.TestCase):
    def test_linear_solve(self):
        ls_entries_orth = [[0, 0, 0, 1],
                           [0, 0, 1, 0],
                           [1, 0, 0, 0],
                           [0, 1, 0, 0]]

        A = OrthogonalMatrix(ls_entries_orth)

        b = Matrix(ls_entries=[[1], [2], [3], [4]])

        x = Matrix(ls_entries=[[3], [4], [2], [1]])
        self.assertEqual(A.solve_linear_system(b), x)
