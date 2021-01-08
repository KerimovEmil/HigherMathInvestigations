from matrix.symmetric_matrix import SymmetricMatrix
from matrix.basic_matrix import MatrixError
import unittest


class HankelMatrix(SymmetricMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)
        if not self.is_hankel():
            raise MatrixError('Not a hankel matrix.')

    def hankel_transform(self):
        """
        Hankel Transform of the Matrix

        Returns:
            <list> sequence of Hankel determinants
        """
        hankel_transform_res = []
        for i in range(1, self.size + 1):
            hankel_transform_res.append(self.matrix_factory(self[:i, : i]).determinant())

        return hankel_transform_res


class TestSquareMatrix(unittest.TestCase):
    def test_hankel_transform(self):
        ls_entries = [[1, 1, 2, 5],
                      [1, 2, 5, 14],
                      [2, 5, 14, 42],
                      [5, 14, 42, 132]]

        A = HankelMatrix(ls_entries)

        self.assertEqual(A.hankel_transform(), [1, 1, 1, 1])
