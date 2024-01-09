from hashlib import sha256
import random

# define OWF {0,1}^n -> {0,1}^n


def o(value: int, f: callable) -> int:
    value_hash = f(str(value).encode('UTF-8')).hexdigest()
    value_hash_int = int(value_hash, 16)
    return value_hash_int


# compute public key using OWF
def generate_public_key(x: list, f: callable, K: int) -> list:
    y = []
    for x_tuple in x:
        y_tuple = (o(x_tuple[0], f),
                   o(x_tuple[1], f))
        y.append(y_tuple)
    return y


# select cryptographic hash function
f = sha256

# set security parameter n to the digest size of the hash function in bits
n = f("".encode("UTF-8")).digest_size * 8

# calculate maximum decimal number that can be represented with n bits
K = pow(2, n) - 1

# calculate private key with values in decimal notation
private = [(random.randint(0, K), random.randint(0, K)) for _ in range(n)]

# calculate public key
public = generate_public_key(private, f, K)

#
# SIGNATURE
#

# get message
M = input("Enter message:")
# hash message
h = f(M.encode('UTF-8')).hexdigest()
# cast hash to binary representation
h_bin = bin(int(h, 16))[2:]

# add missing zeros at the beginning of the binary representation
if int(h[0], 16) < 2:
    h_bin = '000' + h_bin
elif int(h[0], 16) < 4:
    h_bin = '00' + h_bin
elif int(h[0], 16) < 8:
    h_bin = '0' + h_bin

signature = []
# pick values from private key corresponding to h to create signature
for index in range(len(h_bin)):
    signature.append(private[index][int(h_bin[index])])

#
# VERIFY
#

verify = True

# compare hashed values of signature with values of public key
for index in range(len(h_bin)):
    if o(signature[index], f) != public[index][int(h_bin[index])]:
        verify = False
        break

if verify:
    print("Verification successfull.")
else:
    print("Verification failed.")
