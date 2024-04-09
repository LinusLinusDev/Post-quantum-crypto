import numpy as np
import sympy as sp
import random
import math


# print readable values
np.set_printoptions(suppress=True)

seed = 0
random.seed(seed)


class GGH:
    # B = "good" lattice basis as private key
    # n = dimension of the lattice
    # p = max absolute value of components of error
    # H = "bad" lattice basis as public key using the Hermite Normal Form of the private key
    def __init__(self, B, p: int):
        self.linear_independence(B)
        self.__B = B
        self.__n = self.__B.shape[0]
        self.__p = p
        self.__H = self.HNF()[0]

    # check if columnvectors in basis are linearly independent
    def linear_independence(self, matrix):
        _, inds = sp.Matrix(matrix).rref()
        if len(inds) != np.shape(matrix)[0]:
            raise Exception(
                "The columnvectors of the basis B have to be linearly independent.")

    def get_public(self):
        return self.__H

    def get_p(self):
        return self.__p

    # extended euclidian algorithm
    def gcd_ext(self, r0: int, r1: int) -> tuple:
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
    def HNF(self):
        # perform Nemhauser/Wolsey algorithm as described here (https://kola.opus.hbz-nrw.de/frontdoor/deliver/index/docId/211/file/Studienarbeit_Kerstin_Susewind.pdf) on page 34.
        H = self.__B
        n = self.__n
        U = np.identity(n)
        i = 0
        while True:
            # Step 1
            for j in range(i+1, n):
                # Step 2
                if H[i, j] != 0:
                    r, p, q = self.gcd_ext(H[i, i], H[i, j])
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

    def encrypt(self, x, pk, p):
        # map message to latticepoint
        m = pk.dot(x)

        # generate random noise vector with magnitude 2
        e = np.array([random.randint(-p, p)
                     for _ in range(self.__n)])

        # encrypt by adding noise vector to lattice point
        c = m + e

        return c

    def decrypt(self, c):
        # multiply with inverted B
        B_inv_c = np.linalg.inv(self.__B).dot(c)

        # round values to closest integer
        B_inv_c_rounded = B_inv_c.round()

        # multiply with B to recover m
        recovered_m = self.__B.dot(B_inv_c_rounded)

        # solve for H to recover x
        recovered_x = np.linalg.inv(self.__H).dot(recovered_m)

        return recovered_x

    def decrypt_with_H(self, c):
        # multiply with inverted B
        H_inv_c = np.linalg.inv(self.__H).dot(c)

        # round values to closest integer
        H_inv_c_rounded = H_inv_c.round()

        # multiply with B to recover m
        recovered_m = self.__H.dot(H_inv_c_rounded)

        # solve for H to recover x
        recovered_x = np.linalg.inv(self.__H).dot(recovered_m)

        return recovered_x


# since this one is hard-coded, you have to edit it here if you want to try another one
# it has to be quadratic and the columnvectors should be linearly independent
base = np.array([[4, -2, 1, 0],
                 [0, -1, 5, 2],
                 [-1, 6, 1, -1],
                 [0, 1, -1, 6]], )

X = GGH(base, 2)

message = np.array([3, -5, 6, -12])

print(f"Public key:")
print(X.get_public())
print()

print(f"Message: {message}")
print()

encrypted_message = X.encrypt(message, X.get_public(), X.get_p())

print(f"Encrypted message: {encrypted_message}")
print()

recovered_message = X.decrypt(encrypted_message)

print(f"Recovered message: {recovered_message}")
print()

recovered_message_public = X.decrypt_with_H(encrypted_message)

print(
    f"Recovered message using public key instead of private key: {recovered_message_public}")
print()
