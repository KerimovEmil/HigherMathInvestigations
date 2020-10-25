from matrix.square_matrix import SquareMatrix
import numpy as np

ls_entries = [[16, 14, 11, 18, 11],
              [15, 14, 8, 15, 10],
              [15, 10, 4, 7, 15],
              [7, 13, 9, 19, 9],
              [12, 4, 19, 8, 9]]

A = SquareMatrix(ls_entries)
A_numpy = np.array(ls_entries)
# print(A.char_eqn_berkowitz().ls_entries[0])
# print(np.poly(A_numpy))

eigs = A.eigenvalues()
# coeff = [1, -15, -53, 0, 42]
# print(np.roots(coeff))
np_eigs, np_eig_vects = np.linalg.eig(A_numpy)

result = all([abs(i - j) <= 1e-5 for i, j in zip(sorted(list(np_eigs)), sorted(eigs))])
print(result)