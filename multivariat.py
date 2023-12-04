import random
import galois
import numpy as np
import sympy as sp


class UOV:
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

    def get_public(self):
        public = []
        substitutions = {}
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        for x in range(self.__n):
            substitutions[variables[x]] = S[x]
        for m in range(self.__o):
            equation = self.__private[m].subs(substitutions)
            equation = sp.expand(equation)
            equation %= self.__K
            public.append(equation)
        return public

    def generate_quadratic(self):
        equation = [np.random.randint(
            self.__K, size=(self.__n, self.__n))] + [np.random.randint(self.__K, size=self.__n)] + [random.randint(0, self.__K - 1)]

        # lower triangle not required
        for i in range(1, self.__n):
            for j in range(i):
                equation[0][i, j] = 0

        # no summands containing more than one oil-variable as factor
        for i in range(self.__o):
            for j in range(i, self.__o):
                equation[0][i, j] = 0
        return equation

    def generate_linear(self):
        equation = [np.random.randint(
            self.__K, size=self.__n)] + [random.randint(0, self.__K - 1)]
        return equation

    def generate_private(self):
        variables = sp.symbols(f"x':{self.__n}", integer=True)
        private = []
        for _ in range(self.__o):
            equation_matrix = self.generate_quadratic()
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

    def generate_S(self):
        variables = sp.symbols(f"x:{self.__n}", integer=True)
        S = []
        for _ in range(self.__n):
            equation_matrix = self.generate_linear()
            equation = 0
            for i in range(self.__n):
                if equation_matrix[0][i] > 0:
                    equation += equation_matrix[0][i]*variables[i]
            if equation_matrix[1] > 0:
                equation += equation_matrix[1]
            S.append(equation)
        return S

    def sign(self, Y: list):
        GF = galois.GF(self.__K)
        v_variables = sp.symbols(f"x'{self.__o}:{self.__n}", integer=True)
        o_variables = sp.symbols(f"x':{self.__o}", integer=True)
        substitutions = {}

        vinegar = GF([random.randint(0, self.__K - 1)
                     for _ in range(self.__v)])
        for i in range(self.__v):
            substitutions[v_variables[i]] = vinegar[i]
        equations = []

        for m in range(self.__o):
            equation = self.__private[m].subs(substitutions)
            equation -= Y[m]
            equations.append(equation)

        equations_x, equations_y = sp.linear_eq_to_matrix(
            equations, o_variables)

        A = GF.Zeros((self.__o, self.__o))
        b = GF.Zeros(self.__o)
        for i in range(self.__o):
            for j in range(self.__o):
                A[i, j] = int(equations_x[i, j] % self.__K)
            b[i] = int(equations_y[i] % self.__K)

        x = np.linalg.solve(A, b)

        return np.concatenate([x, vinegar])

        """print(temp_y)
       # A = GF(equations_x)
        b = GF(temp_y)"""


X = UOV(2, 3)
print(X.get_private())
print(X.sign([0, 1]))
