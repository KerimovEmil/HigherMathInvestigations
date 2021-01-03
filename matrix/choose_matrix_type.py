from matrix.square_matrix import SquareMatrix
from matrix.symmetric_matrix import SymmetricMatrix
from matrix.hankel_matrix import HankelMatrix
from matrix.basic_matrix import Matrix, MatrixError


# main link: https://www.mins.ee.ethz.ch/teaching/ha/handouts/linalg3p.pdf


class MatrixFactory:  # todo consider making this a function instead of a class
    def __call__(self, ls_entries=None):
        """Returns the most relevant matrix type that exists."""
        m_obj = Matrix(ls_entries)  # todo: consider using Matrix.__subclasses__()

        if m_obj.is_square():
            m_obj = SquareMatrix(ls_entries)
            if m_obj.is_symmetric():
                m_obj = SymmetricMatrix(ls_entries)
                if m_obj.is_hankel():
                    m_obj = HankelMatrix(ls_entries)
        return m_obj
