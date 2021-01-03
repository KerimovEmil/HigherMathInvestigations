from matrix.symmetric_matrix import SymmetricMatrix
from matrix.basic_matrix import MatrixError


class HankelMatrix(SymmetricMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)
        if not self.is_hankel():
            raise MatrixError('Not a hankel matrix.')
