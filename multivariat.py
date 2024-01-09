import random
import galois
import numpy as np
import sympy as sp

# different random maps and Vinegar variables for seed = -1, same random maps and Vinegar variables for seed > -1
seed = -1

if seed >= 0:
    random.seed(seed)
    np.random.seed(seed)


class UOV:
    # o = number of Oil variables
    # v = number of Vinegar variables
    # K = modulus of field
    # n = number of total variables
    # S = affine map
    # private = quadratic map with trapdoor
    def __init__(self, o: int, v: int, K: int = 2):
        self.__o = o
        self.__v = v
        self.__K = K
        self.__n = o + v
        self.__S = self.generate_S()
        self.__private = self.generate_private()

    def get_S(self):
        return self.__S

    def get_private(self):
        return self.__private

    # Get composition of affine and quadratic map as public key
    def get_public(self):
        public = []
        substitutions = {}
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        for x in range(self.__n):
            substitutions[variables[x]] = self.__S[x]
        # substitute each variable in each row of the quadratic map by the corresponding row of the affine map and then simplify the equations
        for m in range(self.__o):
            equation = self.__private[m].subs(substitutions)
            equation = sp.expand(equation, modulus=self.__K)
            public.append(equation)
        return public

    # generate single multivariate quadratic equation with n variables
    def generate_quadratic(self):
        # matrix stands for the quadratic part, list for the linear part, integer for the constant part
        equation = [np.random.randint(
            self.__K, size=(self.__n, self.__n))] + [np.random.randint(self.__K, size=self.__n)] + [random.randint(0, self.__K - 1)]

        # lower triangle not required since multiplication is commutative
        for i in range(1, self.__n):
            for j in range(i):
                equation[0][i, j] = 0

        # no summands containing more than one oil-variable as factor
        for i in range(self.__o):
            for j in range(i, self.__o):
                equation[0][i, j] = 0
        return equation

    # generate affine map with n variables and equations, has to be invertible over finite field K
    def generate_linear(self):
        GF = galois.GF(self.__K)
        while True:
            matrix = GF.Random((self.__n, self.__n))
            # stop generating if invertible map is found
            if np.linalg.det(matrix) != 0:
                break
        # generate constants
        vector = np.random.randint(self.__K, size=self.__n)
        return [matrix, vector]

    # generate o quadratic equations with n variables and create sympy expressions
    def generate_private(self):
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        private = []
        for _ in range(self.__o):
            # generate quadratic equation
            equation_matrix = self.generate_quadratic()
            # transform to sympy expression
            equation = 0
            for i in range(self.__n):
                for j in range(self.__n):
                    if equation_matrix[0][i, j] > 0:
                        equation += equation_matrix[0][i,
                                                       j]*variables[i]*variables[j]
            for i in range(self.__n):
                if equation_matrix[1][i] > 0:
                    equation += equation_matrix[1][i]*variables[i]
            if equation_matrix[2] > 0:
                equation += equation_matrix[2]
            private.append(equation)
        return private

    # generate affine map and create sympy expression
    def generate_S(self):
        variables = sp.symbols(f"x:{self.__n}", integer=True)
        S = []
        # generate linear equations
        equations_matrix = self.generate_linear()
        # transform to sympy expressions
        for i in range(self.__n):
            equation = 0
            for j in range(self.__n):
                if equations_matrix[0][i, j] > 0:
                    equation += int(equations_matrix[0][i, j])*variables[j]
            if equations_matrix[1][i] > 0:
                equation += equations_matrix[1][i]
            S.append(equation)
        return S

    # sign list of o integers mod K
    def sign(self, Y: list):
        GF = galois.GF(self.__K)
        A = GF.Zeros((self.__o, self.__o))
        b = GF.Zeros(self.__o)

        v_variables = sp.symbols(f"x'{self.__o}:{self.__n}", integer=True)
        o_variables = sp.symbols(f"x':{self.__o}", integer=True)
        variables = sp.symbols(f"x:{self.__n}", integer=True)

        substitutions = {}
        failures = 0

        # inverting quadratic map
        while True:
            # set vinegar variables randomly and substitute them in private system
            vinegar = GF([random.randint(0, self.__K - 1)
                          for _ in range(self.__v)])
            for i in range(self.__v):
                substitutions[v_variables[i]] = vinegar[i]
            equations = []

            for m in range(self.__o):
                equation = self.__private[m].subs(substitutions)
                equation = sp.expand(equation, modulus=self.__K)
                equations.append(equation)

            # cast linear system to GF matrices to solve them over finite field K
            equations_x, equations_y = sp.linear_eq_to_matrix(
                equations, o_variables)

            for i in range(self.__o):
                for j in range(self.__o):
                    A[i, j] = int(equations_x[i, j])
                # add values from Y
                b[i] = (int(equations_y[i]) + Y[i]) % self.__K

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
        A = GF.Zeros((self.__n, self.__n))
        b = GF.Zeros(self.__n)
        equations_x, equations_y = sp.linear_eq_to_matrix(
            self.__S, variables)
        for i in range(self.__n):
            for j in range(self.__n):
                A[i, j] = int(equations_x[i, j])
            b[i] = (int(equations_y[i]) + int(y_new[i])) % self.__K

        return np.linalg.solve(A, b)


X = UOV(2, 2)
document = [1, 0]

print(f"Private system: {X.get_private()}")
print()
print(f"Private map S: {X.get_S()}")
print()
print(f"Public system: {X.get_public()}")
print()
print(f"Document: {document}")
print(f"Signature: {X.sign(document)}")
