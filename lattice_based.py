import numpy as np
import sympy as sp
import random
import math

# print readable values
np.set_printoptions(suppress=True)

# different random errors for seed = -1, same random errors for seed > -1
seed = -1

if seed >= 0:
    random.seed(seed)


# check if columnvectors in basis are linear independent
def linear_independence(matrix):
    _, inds = sp.Matrix(matrix).rref()
    if len(inds) == np.shape(matrix)[0]:
        return 1
    else:
        return 0


# extended euclidian algorithm
def gcd_ext(r0: int, r1: int) -> tuple:
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


# lower triangular basis using Nemhauser/Wolsey algorithm
def HNF(A):
    # perform Nemhauser/Wolsey algorithm as described here (https://kola.opus.hbz-nrw.de/frontdoor/deliver/index/docId/211/file/Studienarbeit_Kerstin_Susewind.pdf) on page 34.
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


# generate noise-vector
def get_errors(A):
    magnitudes = []
    # compute magnitudes of all columnvectors
    for vector in A.T:
        magnitudes.append(np.linalg.norm(vector))
    # generate random noise-vector
    error = np.array([random.uniform(-1, 1) for _ in range(np.shape(A)[0])])
    # scale noise-vector to 30%-50% magnitude of shortest columnvector in basis
    error *= min(magnitudes) * random.uniform(0.3, 0.5) / \
        np.linalg.norm(error)
    return error


# Babai's rounding procedure
def recover_v(c, B):
    # multiply with inverted B
    B_inv_c = np.linalg.inv(B).dot(c)

    # round values to closest integer
    B_inv_c_rounded = B_inv_c.round()

    # multiply with B to recover v
    recovered_v = B.dot(B_inv_c_rounded)

    return recovered_v


# since this one is hard-coded, you have to edit it here if you want to try another one
# it has to be quadratic and the columnvectors should be linear independent
B = np.array([[4, -2, 1, 0],
              [0, -1, 5, 2],
              [-1, 6, 1, -1],
              [0, 1, -1, 6]])


class GHH:
    # B = "good" lattice basis as private key
    def __init__(self, B):
        self.__B = B
        self.__n = self.__B.shape[0]


if linear_independence(B) == 0:
    print("The columnvectors of the basis B are not linear independent.")
else:
    # "bad" lattice basis as public key using the Hermite Normal Form of the public key
    H, U = HNF(B)

    # get dimension of basis
    dim = np.shape(B)[0]

    # get message x and compute lattice point m
    message = input(
        f"Type in {dim} integers es message, seperated by a comma:")
    message = message.split(",")
    x = [int(x) for x in message]
    m = H.dot(np.array(x))

    # short noise vector
    e = get_errors(B)

    # encrypted lattice point by adding noise
    c = m + e

    print()
    print("B = ")
    print(B)
    print("H = ")
    print(H)
    print()
    print(f"x = {x}")
    print(f"m = {m}")
    print(f"e = {e}")
    print(f"c = {c}")

    # recover lattice point and message using good and bad base
    recovered_B = recover_v(c, B)
    recovered_H = recover_v(c, H)
    print(f"recovered point with good base B = {recovered_B}")
    print(
        f"recovered message with good base B = {np.linalg.inv(H).dot(recovered_B)}")
    print(f"recovered point with bad base H = {recovered_H}")
    print(
        f"recovered message with bad base H = {np.linalg.inv(H).dot(recovered_H)}")
