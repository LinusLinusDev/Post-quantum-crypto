import random
import numpy as np
from functions import HNF, recover_v


# "good" lattice basis as private key
B = np.array([[1, 2, 3, 4, 5],
              [2, 3, 4, 5, 1],
              [3, 4, 5, 1, 2],
              [4, 5, 1, 2, 3],
              [5, 1, 2, 3, 4]])

# "bad" lattice basis as public key using the Hermite Normal Form of the public key
H, U = HNF(B)

# get message x and compute lattice point m
message = input("Type in 5 integers es message, seperated by a comma:")
message = message.split(",")
x = [int(x) for x in message]
m = H.dot(np.array(x))

# short noise vector
errors = [random.uniform(-1.5, 1.5) for _ in range(5)]
e = np.array(errors)

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
print(f"recovered message with good base B = {np.linalg.inv(H).dot(recovered_B)}")
print(f"recovered point with bad base H = {recovered_H}")
print(f"recovered message with bad base H = {np.linalg.inv(H).dot(recovered_H)}")
