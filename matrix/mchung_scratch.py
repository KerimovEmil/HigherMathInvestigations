from matrix.square_matrix import SquareMatrix
import numpy as np

ls_entries = [[1, 2, 1, 5], [1, 3, 2, 5], [1, 5, 2, 5], [6, 7, 8, 9]]

A = SquareMatrix(ls_entries=ls_entries)
A_numpy = np.array(ls_entries)
# print(A.char_eqn_berkowitz().ls_entries[0])
# print(np.poly(A_numpy))

print(A.eigenvalues())
# coeff = [1, -15, -53, 0, 42]
# print(np.roots(coeff))
print(np.linalg.eig(A_numpy))