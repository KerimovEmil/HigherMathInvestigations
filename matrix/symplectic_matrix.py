from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError
import unittest


class SymplecticMatrix(SquareMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

        if not self.is_symplectic():
            raise MatrixError('Not a symplectic matrix.')

    def determinant(self):
        # determinant of a symplectic matrix is 1
        return 1

    def inverse(self):
        """ Get the inverse of the symplectic matrix """
        # self is a 2n x 2n matrix
        # Use the typical block matrix:
        # [[0  I_n],
        # [-I_n 0]]
        n = self.size // 2  # be an int
        omega = self.block_matrix(top_left=self.zero_matrix(n, n),
                                  top_right=self.identity(n),
                                  bottom_left=-self.identity(n),
                                  bottom_right=self.zero_matrix(n, n))

        # use parent method since omega is
        return SquareMatrix.inverse(omega) * self.transpose() * omega


class TestSymplecticMatrix(unittest.TestCase):
    def test_inverse(self):
        import numpy as np

        A = SymplecticMatrix(ls_entries=[[1, 0, 0, 1],
                                         [0, 1, 1, 0],
                                         [0, 0, 1, 0],
                                         [0, 0, 0, 1]])

        B = SymplecticMatrix(ls_entries=[[0, 1, 0, 1],
                                         [1, 0, 1, 0],
                                         [0, 0, 0, 1],
                                         [0, 0, 1, 0]])

        self.assertTrue(np.allclose(A.inverse().ls_entries, np.linalg.inv(A.ls_entries)))
        self.assertTrue(np.allclose(B.inverse().ls_entries, np.linalg.inv(B.ls_entries)))
