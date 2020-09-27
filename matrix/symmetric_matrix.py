from matrix.square_matrix import SquareMatrix


class SymmetricMatrix(SquareMatrix):
    def transpose(self):
        return SymmetricMatrix(ls_entries=self.ls_entries)
