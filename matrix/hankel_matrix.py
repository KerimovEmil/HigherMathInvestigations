from matrix.symmetric_matrix import SymmetricMatrix


class HankelMatrix(SymmetricMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)

