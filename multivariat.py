import random
import galois
import numpy as np
import sympy as sp

seed = 6
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
    def __init__(self, o: int, v: int, p: int):
        self.__o = o
        self.__v = v
        self.__n = self.__o + self.__v
        self.__p = p
        self.__K = galois.GF(self.__p)
        self.__S = self.generate_S()
        self.__private = self.generate_private()
        self.__public = self.generate_public()

    def get_public(self):
        return self.__public

    def get_private(self):
        return self.__private

    def get_S(self):
        return self.__S

    def to_sympy(self, system, quadratic: bool):
        variables = sp.symbols(f"x:{self.__n}", integer=True)
        sym = []
        if quadratic:
            for o in range(self.__o):
                function = 0
                for i in range(self.__n):
                    for j in range(self.__n):
                        function += int(system[o][0]
                                        [i, j])*variables[i]*variables[j]
                for i in range(self.__n):
                    function += int(system[o][1][i])*variables[i]
                function += system[o][2]
                sym.append(function)
        else:
            for i in range(self.__n):
                function = 0
                for j in range(self.__n):
                    function += int(system[0][i, j])*variables[j]
                function += system[1][i]
                sym.append(function)
        return sym

    # Get composition of affine and quadratic map as public key
    def generate_public(self):
        D = self.__S[0]
        e = self.__S[1]
        public = []

        for i in range(self.__o):
            # calculate composition
            A = self.__private[i][0]
            b = self.__private[i][1]
            c = self.__private[i][2]
            function = [D.transpose().dot(A).dot(D),
                        (e.transpose().dot(A.transpose()+A)+b.transpose()).dot(D),
                        (e.transpose().dot(A)+b.transpose()).dot(e)+c]

            # summarize quadratic coefficients in upper triangle
            for j in range(self.__n):
                for k in range(j + 1, self.__n):
                    function[0][j, k] = function[0][j, k] + function[0][k, j]
                    function[0][k, j] = 0

            public.append(function)

        return public

    # generate affine map with n variables and equations, has to be invertible over finite field K
    def generate_S(self):
        # generate invertible matrix
        D = self.__K.Random((self.__n, self.__n),
                            seed=random.randint(0, 1000))
        # generate constants
        e = self.__K.Random(self.__n, seed=random.randint(0, 1000))

        while np.linalg.det(D) == 0:
            D = self.__K.Random(
                (self.__n, self.__n), seed=random.randint(0, 1000))

        return [D, e]

    # generate o quadratic equations with n variables
    def generate_private(self):
        private = []
        for _ in range(self.__o):
            # matrix stands for the quadratic part, list for the linear part, integer for the constant part
            function = [self.__K.Random((self.__n, self.__n), seed=random.randint(0, 1000))] + [self.__K.Random(
                self.__n, seed=random.randint(0, 1000))] + [self.__K.Random(seed=random.randint(0, 1000))]

            # set upper left o x o Part of matrix to 0, since products with more than one oil-factor are not allowed
            for i in range(self.__n):
                for j in range(max(self.__o, i)):
                    function[0][i, j] = 0
            private.append(function)
        return private

    # sign list of o integers mod K
    def sign(self, y: list):
        failures = 0

        # inverting quadratic map
        while True:
            # set vinegar variables randomly and substitute them in private system
            x_v = self.__K.Random(self.__v, seed=random.randint(0, 1000))

            linear = []
            constant = []

            for i in range(self.__o):
                M = self.__private[i][0][np.ix_(
                    range(self.__o), range(self.__o, self.__n))]
                N = self.__private[i][0][np.ix_(
                    range(self.__o, self.__n), range(self.__o, self.__n))]
                m = self.__private[i][1][np.ix_(range(self.__o))]
                n = self.__private[i][1][np.ix_(range(self.__o, self.__n))]
                c = self.__private[i][2]

                linear.append(x_v.transpose().dot(M.transpose())+m.transpose())
                constant.append(
                    (x_v.transpose().dot(N)+n.transpose()).dot(x_v)+c)

            linear_all = self.__K(linear)
            constant_all = self.__K(constant)

            # try to solve system, if not possible set different vinegar variables
            try:
                x_o = np.linalg.solve(linear_all, self.__K(y) - constant_all)
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
        x_line = np.concatenate([x_o, x_v])

        # invert affine map over finite field
        return np.linalg.solve(self.__S[0], x_line - self.__S[1])

    # verify signature
    def verify(self, x, y, P):
        result = []

        for i in range(len(P)):
            A = P[i][0]
            b = P[i][1]
            c = P[i][2]
            result.append(x.transpose().dot(A).dot(
                x)+b.transpose().dot(x)+c)

        return result == y


X = UOV(4, 4, 2)

print()
print(f"Public system: {X.to_sympy(X.get_public(),True)}")
print()
print(f"Private systen: {X.to_sympy(X.get_private(),True)}")
print()
print(f"S: {X.to_sympy(X.get_S(),False)}")
print()

document = [1, 0, 1, 1]

print(f"Document: {document}")
print()

signature = X.sign(document)

print(f"Signature: {signature}")
print()

if X.verify(signature, document, X.get_public()):
    print("The verification was successful.")
else:
    print("The verification was not successful.")
