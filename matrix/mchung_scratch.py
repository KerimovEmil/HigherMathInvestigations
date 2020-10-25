from matrix.square_matrix import SquareMatrix
import numpy as np

ls_entries = [[1, 2, 1, 5], [1, 3, 2, 5], [1, 5, 2, 5], [6, 7, 8, 9]]

A = SquareMatrix(ls_entries=ls_entries)
A_numpy = np.array(ls_entries)
print(A.char_eqn_berkowitz())
print(np.poly(A_numpy))
