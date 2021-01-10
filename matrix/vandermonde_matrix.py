from matrix.basic_matrix import MatrixError, Matrix
from matrix.square_matrix import SquareMatrix

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