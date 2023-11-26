import numpy as np
from lattice_based_functions import linear_independence, HNF, get_errors, recover_v

np.set_printoptions(suppress=True)

# "good" lattice basis as private key

# since this one is hard-coded, you have to edit it here if you want to try another one
# it has to be quadratic and the columnvectors should be linear independent
B = np.array([[1, 2, 3, 4, 5, 6],
              [2, 3, 4, 5, 6, 1],
              [3, 4, 5, 6, 1, 2],
              [4, 5, 6, 1, 2, 3],
              [5, 6, 1, 2, 3, 4],
              [6, 1, 2, 3, 4, 5]])

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
