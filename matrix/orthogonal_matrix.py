from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError
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
            b: constant matrix

        Returns:
            x: vector solution of linear system
        """
        return self.transpose().__mul__(b)

class TestOrthogonalMatrix(unittest.TestCase):
    def test_linear_solve(self):
