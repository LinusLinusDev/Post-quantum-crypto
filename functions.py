import math
import numpy as np


def linear_independence(matrix):
    independence = 1
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i != j:
                inner_product = np.inner(
                    matrix[:, i],
                    matrix[:, j]
                )
                norm_i = np.linalg.norm(matrix[:, i])
                norm_j = np.linalg.norm(matrix[:, j])

                if np.abs(inner_product - norm_j * norm_i) < 1E-5:
                    independence = 0
    return independence


def gcd_ext(r0: int, r1: int) -> tuple:  # extended euclidian algorithm
    if r0 < 0:
        vz0 = -1
        r0 = -r0
    else:
        vz0 = 1
    if r1 < 0:
        vz1 = -1
        r1 = -r1
    else:
        vz1 = 1

    x0, x1, y0, y1 = 1, 0, 0, 1

    while True:
        if r1 == 0:
            return r0, vz0 * x0, vz1 * y0
        q = r0 // r1
        r0, r1 = r1, r0 % r1
        x0, x1 = x1, x0 - x1 * q
        y0, y1 = y1, y0 - y1 * q


def recover_v(c, B):  # Babai's rounding procedure
    # multiply with inverted B
    B_inv_c = np.linalg.inv(B).dot(c)

    # round values to closest integer
    B_inv_c_rounded = B_inv_c.round()

    # multiply with B to recover v
    recovered_v = B.dot(B_inv_c_rounded)

    return recovered_v


def HNF(A):  # lower triangular basis using Nemhauser/Wolsey algorithm
    H = A
    n = H.shape[0]
    U = np.identity(n)
    i = 0
    while True:
        # Step 1
        j = i + 1
        while j <= n - 1:
            # Step 2
            if H[i, j] != 0:
                r, p, q = gcd_ext(H[i, i], H[i, j])
                temp = np.identity(n)
                temp[i, i] = p
                temp[j, i] = q
                temp[i, j] = -H[i, j] / r
                temp[j, j] = H[i, i] / r
                H = H.dot(temp)
                U = U.dot(temp)
            j = j + 1
        # Step 3
        if H[i, i] < 0:
            temp = np.identity(n)
            temp[i, i] = -1
            H = H.dot(temp)
            U = U.dot(temp)
        j = 0
        while True:
            # Step 4
            if j == i:
                if i == n - 1:
                    return (H, U)
                else:
                    i = i + 1
                    break
            temp = np.identity(n)
            temp[i, j] = -math.ceil(H[i, j] / H[i, i])
            H = H.dot(temp)
            U = U.dot(temp)
            if j == i - 1:
                i = i + 1
                if i > n - 1:
                    return (H, U)
                else:
                    break
            if j < i - 1:
                j = j + 1
