from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError


class SkewSymmetric(SquareMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

        if not self.is_skew_symmetric():
            raise MatrixError('Not a skew symmetric matrix.')

    def transpose(self):
        return -self

