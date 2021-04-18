from matrix.square_matrix import SquareMatrix
from matrix.symmetric_matrix import SymmetricMatrix
from matrix.hankel_matrix import HankelMatrix
from matrix.basic_matrix import Matrix, MatrixError
from matrix.vandermonde_matrix import VandermondeMatrix, SquareVandermondeMatrix


# main link: https://www.mins.ee.ethz.ch/teaching/ha/handouts/linalg3p.pdf


class MatrixFactory:  # todo consider making this a function instead of a class
    def get_matrix_class(self, m_obj, ls_entries):
        """ recursively loop through the subclasses to get the most relevant type"""
        first_class = m_obj.__class__
        subclass_ls = first_class.__subclasses__()
        if subclass_ls:
            for subcls in subclass_ls:
                try:
                    m_obj = subcls(ls_entries)
                except MatrixError:
                    pass

        if first_class == m_obj.__class__:
            return m_obj
        else:
            return self.get_matrix_class(m_obj, ls_entries)

    def __call__(self, ls_entries=None):
        """Returns the most relevant matrix type that exists."""
        m_obj = Matrix(ls_entries)
        m_obj = self.get_matrix_class(m_obj, ls_entries)

        return m_obj
