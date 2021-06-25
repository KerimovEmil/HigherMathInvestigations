from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError


class SymplecticMatrix(SquareMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

        if not self.is_symplectic():
            raise MatrixError('Not a symplectic matrix.')

