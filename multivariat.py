import random
import galois
import numpy as np
import sympy as sp

# different random maps and Vinegar variables for seed = -1, same random maps and Vinegar variables for seed > -1
seed = 6

if seed >= 0:
    random.seed(seed)


class UOV:
    # o = number of Oil variables
    # v = number of Vinegar variables
    # n = number of total variables
    # m = modulus of finite field
    # K = finite field
    # S = affine map
    # private = quadratic map with trapdoor
    # public = composition of affine and quadratic map
    def __init__(self, o: int, v: int, m: int):
        self.__o = o
        self.__v = v
        self.__n = self.__o + self.__v
        self.__m = m
        self.__K = galois.GF(self.__m)
        self.__S = self.generate_S()
        self.__private = self.generate_private()
        self.__public = self.generate_public()

    def get_public(self):
        return self.__public

    def get_private(self):
        return self.__private

    def get_S(self):
        return self.__S

    def to_sympy(self, A, quadratic: bool):
        variables = sp.symbols(f"x:{self.__n}", integer=True)
        sym = []
        if quadratic:
            for o in range(self.__o):
                equation = 0
                for i in range(self.__n):
                    for j in range(self.__n):
                        equation += int(A[o][0]
                                        [i, j])*variables[i]*variables[j]
                for i in range(self.__n):
                    equation += int(A[o][1][i])*variables[i]
                equation += A[o][2]
                sym.append(equation)
        else:
            for i in range(self.__n):
                equation = 0
                for j in range(self.__n):
                    equation += int(A[0][i, j])*variables[j]
                equation += A[1][i]
                sym.append(equation)
        return sym

    # Get composition of affine and quadratic map as public key
    def generate_public(self):
        S = self.__S[0]
        r = self.__S[1]
        public = []

        for i in range(self.__o):
            # calculate composition
            A = self.__private[i][0]
            b = self.__private[i][1]
            c = self.__private[i][2]
            equation = [S.transpose().dot(A).dot(S), (r.transpose().dot(A.transpose(
            )+A)+b.transpose()).dot(S), (r.transpose().dot(A)+b.transpose()).dot(r)+c]
            public.append(equation)

            # summarize quadratic coefficients in upper triangle
            for j in range(self.__n):
                for k in range(j + 1, self.__n):
                    equation[0][j, k] = equation[0][j, k] + equation[0][k, j]
                    equation[0][k, j] = 0

        return public

    # generate affine map with n variables and equations, has to be invertible over finite field K
    def generate_S(self):
        # generate invertible matrix
        matrix = self.__K.Random((self.__n, self.__n),
                                 seed=random.randint(0, 1000))
        while np.linalg.det(matrix) == 0:
            matrix = self.__K.Random(
                (self.__n, self.__n), seed=random.randint(0, 1000))
        # generate constants
        vector = self.__K.Random(self.__n, seed=random.randint(0, 1000))
        return [matrix, vector]

    # generate o quadratic equations with n variables
    def generate_private(self):
        private = []
        for _ in range(self.__o):
            # matrix stands for the quadratic part, list for the linear part, integer for the constant part
            equation = [self.__K.Random((self.__n, self.__n), seed=random.randint(0, 1000))] + [self.__K.Random(
                self.__n, seed=random.randint(0, 1000))] + [self.__K.Random(seed=random.randint(0, 1000))]

            # set upper left o x o Part of matrix to 0, since products with more than one oil-factor are not allowed
            for i in range(self.__n):
                for j in range(max(self.__o, i)):
                    equation[0][i, j] = 0
            private.append(equation)
        return private

    # sign list of o integers mod K
    def sign(self, Y: list):
        failures = 0

        # inverting quadratic map
        while True:
            # set vinegar variables randomly and substitute them in private system
            v = self.__K.Random(self.__v, seed=random.randint(0, 1000))

            linear = []
            constant = []

            for i in range(self.__o):
                A_1 = self.__private[i][0][np.ix_(
                    range(self.__o), range(self.__o, self.__n))]
                A_2 = self.__private[i][0][np.ix_(
                    range(self.__o, self.__n), range(self.__o, self.__n))]
                b_1 = self.__private[i][1][np.ix_(range(self.__o))]
                b_2 = self.__private[i][1][np.ix_(range(self.__o, self.__n))]
                c = self.__private[i][2]

                linear.append(v.transpose().dot(
                    A_1.transpose())+b_1.transpose())
                constant.append(
                    (v.transpose().dot(A_2)+b_2.transpose()).dot(v)+c)

            linear_all = self.__K(linear)
            constant_all = self.__K(constant)

            # try to solve system, if not possible set different vinegar variables
            try:
                x = np.linalg.solve(linear_all, self.__K(Y) - constant_all)
                break
            except:
                failures += 1
                print(f"Failure {failures}")
                # stop trying after 10 attempts
                if failures == 10:
                    print("Can not sign document.")
                    return -1
                continue

        # safe results and vinegar variables as total result of inverted quadratic system
        y_new = np.concatenate([x, v])

        # invert affine map over finite field
        return np.linalg.solve(self.__S[0], y_new - self.__S[1])

    # verify signature
    def verify(self, signature, M, pk):
        message = []

        for i in range(len(pk)):
            A = pk[i][0]
            b = pk[i][1]
            c = pk[i][2]
            message.append(signature.transpose().dot(A).dot(
                signature)+b.transpose().dot(signature)+c)

        return message == M


X = UOV(8, 8, 19)

print()
print(f"Public system: {X.to_sympy(X.get_public(),True)}")
print()
print(f"Private systen: {X.to_sympy(X.get_private(),True)}")
print()
print(f"S: {X.to_sympy(X.get_S(),False)}")
print()

document = [1, 0, 1, 1, 1, 1, 1, 1]

print(f"Document: {document}")
print()

signature = X.sign(document)

print(f"Signature: {signature}")
print()

if X.verify(signature, document, X.get_public()):
    print("The verification was successful.")
else:
    print("The verification was not successful.")
