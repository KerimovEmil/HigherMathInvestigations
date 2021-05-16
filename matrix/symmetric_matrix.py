from matrix.square_matrix import SquareMatrix
from matrix.basic_matrix import MatrixError
import pandas as pd
import unittest


class SymmetricMatrix(SquareMatrix):
    def __init__(self, ls_entries=None, matrix=None):
        super().__init__(ls_entries, matrix)
        if not self.is_symmetric():
            raise MatrixError('Not a symmetric matrix.')

    def transpose(self):
        return self

    def get_definiteness(self):
        """
        Returns: <str> the definiteness type, i.e. positive_definite, positive_semi_definite, negative_definite,
        negative_semi_definite, or indefinite
        """
        eigenvalues = self.eigenvalues()
        definiteness_dict = {'positive_definite': all(i > 0 for i in eigenvalues),
                             'positive_semi_definite': all(i >= 0 for i in eigenvalues),
                             'negative_definite': all(i < 0 for i in eigenvalues),
                             'negative_semi_definite': all(i <= 0 for i in eigenvalues)}

        # get the true conditions
        df = pd.DataFrame(
            {'definiteness_type': list(definiteness_dict.keys()), 'condition_res': list(definiteness_dict.values())})

        possible_types = df[df['condition_res'] == True]['definiteness_type'].tolist()

        if possible_types:
            if len(
                    possible_types) == 2:  # both the definite and semi-definite case have been flagged as true, then it is definite
                res = 'positive_definite' if 'positive_definite' in possible_types else 'negative_definite'
                return res
            return possible_types[0]

        return 'indefinite'


class TestSymmetricMatrix(unittest.TestCase):
    def test_definiteness(self):
        pos_def = SymmetricMatrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
        pos_semi_def = SymmetricMatrix([[2, 6], [6, 18]])
        in_def = SymmetricMatrix([[1, 3 / 4, 0], [3 / 4, 1, 3 / 4], [0, 3 / 4, 1]])
        neg_def = SymmetricMatrix([[-3, 0, 0], [0, -2, 0], [0, 0, -1]])
        neg_semi_def = SymmetricMatrix([[0, 0], [0, -1]])

        self.assertEqual(pos_def.get_definiteness(), 'positive_definite')
        self.assertEqual(pos_semi_def.get_definiteness(), 'positive_semi_definite')
        self.assertEqual(in_def.get_definiteness(), 'indefinite')
        self.assertEqual(neg_def.get_definiteness(), 'negative_definite')
        self.assertEqual(neg_semi_def.get_definiteness(), 'negative_semi_definite')
