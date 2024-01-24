import random
import galois
import numpy as np
import sympy as sp

# different random maps and Vinegar variables for seed = -1, same random maps and Vinegar variables for seed > -1
seed = -1

if seed >= 0:
    random.seed(seed)


class UOV:
    # o = number of Oil variables
    # v = number of Vinegar variables
    # K = modulus of field
    # n = number of total variables
    # S = affine map
    # private = quadratic map with trapdoor
    def __init__(self, o: int, v: int, m: int = 2):
        self.__o = o
        self.__v = v
        self.__n = o + v
        self.__m = m
        self.__K = galois.GF(m)
        self.__Snum = self.generate_Snum()
        self.__Ssym = self.generate_Ssym()
        self.__private = self.generate_private()
        self.__public = self.generate_public()

    def get_S(self):
        return self.__Ssym

    def get_private(self):
        return self.__private

    def get_public(self):
        return self.__public

    # Get composition of affine and quadratic map as public key
    def generate_public(self):
        public = []
        substitutions = {}
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        for x in range(self.__n):
            substitutions[variables[x]] = self.__Ssym[x]
        # substitute each variable in each row of the quadratic map by the corresponding row of the affine map and then simplify the equations
        for m in range(self.__o):
            equation = self.__private[m].subs(substitutions)
            equation = sp.expand(equation, modulus=self.__m)
            public.append(equation)
        return public

    # generate affine map with n variables and equations, has to be invertible over finite field K
    def generate_Snum(self):
        # generate invertible matrix
        matrix = self.__K.Random((self.__n, self.__n),
                                 seed=random.randint(0, 1000))
        while np.linalg.det(matrix) == 0:
            matrix = self.__K.Random(
                (self.__n, self.__n), seed=random.randint(0, 1000))
        # generate constants
        vector = self.__K.Random(self.__n, seed=random.randint(0, 1000))
        return [matrix, vector]

    # generate o quadratic equations with n variables and create sympy expressions
    def generate_private(self):
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        private = []
        for _ in range(self.__o):
            # matrix stands for the quadratic part, list for the linear part, integer for the constant part
            equation_matrix = [self.__K.Random((self.__n, self.__n), seed=random.randint(0, 1000))] + [self.__K.Random(
                self.__n, seed=random.randint(0, 1000))] + [self.__K.Random(seed=random.randint(0, 1000))]
            # transform to sympy expression
            equation = 0
            # skip lower triangle of matrix since its redundant
            # also skip first o x o Part of matrix, since products with more than one oil-factor are not allowed
            for i in range(self.__n):
                for j in range(max(self.__o, i), self.__n):
                    equation += int(equation_matrix[0]
                                    [i, j])*variables[i]*variables[j]
            for i in range(self.__n):
                equation += int(equation_matrix[1][i])*variables[i]
            equation += equation_matrix[2]
            private.append(equation)
        return private

    # create sympy expression using affine Map Snum
    def generate_Ssym(self):
        variables = sp.symbols(f"x:{self.__n}", integer=True)
        S = []
        # transform to sympy expressions
        for i in range(self.__n):
            equation = 0
            for j in range(self.__n):
                equation += int(self.__Snum[0][i, j])*variables[j]
            equation += self.__Snum[1][i]
            S.append(equation)
        return S

    # sign list of o integers mod K
    def sign(self, Y: list):
        A = self.__K.Zeros((self.__o, self.__o))
        b = self.__K.Zeros(self.__o)

        v_variables = sp.symbols(f"x'{self.__o}:{self.__n}", integer=True)
        o_variables = sp.symbols(f"x':{self.__o}", integer=True)

        substitutions = {}
        failures = 0

        # inverting quadratic map
        while True:
            # set vinegar variables randomly and substitute them in private system
            vinegar = self.__K.Random(self.__v, seed=random.randint(0, 1000))
            for i in range(self.__v):
                substitutions[v_variables[i]] = vinegar[i]
            equations = []

            for m in range(self.__o):
                equation = self.__private[m].subs(substitutions)
                equation = sp.expand(equation, modulus=self.__m)
                equations.append(equation)

            # cast linear system to GF matrices to solve them over finite field K
            equations_x, equations_y = sp.linear_eq_to_matrix(
                equations, o_variables)

            for i in range(self.__o):
                for j in range(self.__o):
                    A[i, j] = int(equations_x[i, j])
                # add values from Y
                b[i] = (int(equations_y[i]) + Y[i]) % self.__m

            # try to solve system, if not possible set different vinegar variables
            try:
                x = np.linalg.solve(A, b)
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
        y_new = np.concatenate([x, vinegar])

        # invert affine map over finite field
        return np.linalg.solve(self.__Snum[0], y_new - self.__Snum[1])

    # verify signature
    def verify(self, X, Y):
        message = []
        substitutions = {}
        variables = sp.symbols(f"x:{self.__n}", integer=True)

        # substitute variables with values of the signature and simplify the equations over finite field
        for x in range(self.__n):
            substitutions[variables[x]] = X[x]
        for m in range(self.__o):
            equation = self.__public[m].subs(substitutions)
            equation = sp.expand(equation, modulus=self.__m)
            message.append(equation)
        return message == Y


X = UOV(4, 4)
document = [1, 0, 0, 1]

print(f"Private system: {X.get_private()}")
print()
print(f"Private map S: {X.get_S()}")
print()
print(f"Public system: {X.get_public()}")
print()
print(f"Document: {document}")

signature = X.sign(document)

print(f"Signature: {signature}")

if X.verify(signature, document):
    print("The verification was successful.")
else:
    print("The verification was not successful.")
