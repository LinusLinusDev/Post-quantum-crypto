from hashlib import sha256
import random


def o(value: int, f: callable, K) -> int:
    value_hash = f(str(value).encode('UTF-8')).hexdigest()
    value_hash_int = int(value_hash, 16)
    return value_hash_int


def generate_public_key(x: list, f: callable, K: int) -> list:
    y = []
    for x_tuple in x:
        y_tuple = (o(x_tuple[0], f, K),
                   o(x_tuple[1], f, K))
        y.append(y_tuple)
    return y


f = sha256
n = f("".encode("UTF-8")).digest_size * 8
K = pow(2, n)

M = input("Enter message:")
h = f(M.encode('UTF-8')).hexdigest()
h_bin = bin(int(h, 16))[2:]
if int(h[0], 16) < 4:
    h_bin = '00' + h_bin
elif int(h[0], 16) < 8:
    h_bin = '0' + h_bin

private = [(random.randint(0, K-1), random.randint(0, K-1)) for _ in range(n)]

public = generate_public_key(private, f, K)

signature = []
for index in range(len(h_bin)):
    signature.append(private[index][int(h_bin[index])])

verify = True

for index in range(len(h_bin)):
    if o(signature[index], f, K) != public[index][int(h_bin[index])]:
        verify = False
        break

if verify:
    print("Verification successfull.")
else:
    print("Verification failed.")
