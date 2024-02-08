from hashlib import sha256
import random

seed = 0
random.seed(seed)


class Lamport:
    # hash = cryptographic hash function used for g and f
    # n = security parameter set to the digest size of the hash function in bits
    # K = maximum decimal number that can be represented with n bits
    # private = private key
    # public = public key
    def __init__(self, hash, n):
        self.__hash = hash
        self.__n = n
        self.__K = pow(2, self.__n) - 1
        self.__private = self.generate_private_key()
        self.__public = self.generate_public_key()

    # setup cryptographic hash function {0,1}^* -> {0,1}^n
    def f(self, message: str) -> str:
        # hash message
        h = self.__hash(message.encode('UTF-8')).hexdigest()
        # cast hash to binary representation
        h_bin = bin(int(h, 16))[2:]
        # add missing zeros at the beginning of the binary representation
        length = self.__hash("".encode("UTF-8")).digest_size * 8
        h_bin = h_bin.zfill(length)

        return h_bin[:self.__n]

    # setup OWF function
    def g(self, value: int) -> int:
        return int(self.f(str(value)), 2)

    # generate private key
    def generate_private_key(self):
        return [(random.randint(0, self.__K), random.randint(0, self.__K)) for _ in range(self.__n)]

    # compute public key using OWF
    def generate_public_key(self):
        y = []
        for x_tuple in self.__private:
            y_tuple = (self.g(x_tuple[0]), self.g(x_tuple[1]))
            y.append(y_tuple)
        return y

    def get_private(self):
        return self.__private

    def get_public(self):
        return self.__public

    def sign(self, y):
        M_bin = self.f(y)

        signature = []
        # pick values from private key corresponding to h to create signature
        for index in range(len(M_bin)):
            signature.append(self.__private[index][int(M_bin[index])])

        return signature

    def verify(self, x, y, pk):
        verify = True
        M_bin = self.f(y)
        V = []
        for s in x:
            V.append(self.g(s))

        # compare hashed values of signature with values of public key
        for index in range(len(M_bin)):
            if V[index] != pk[index][int(M_bin[index])]:
                verify = False
                break

        return verify


X = Lamport(sha256, 16)
document = "Hello World!"

print(f"Private key: {X.get_private()}")
print()
print(f"Public key: {X.get_public()}")
print()
print(f"Document: {document}")
print()

signature = X.sign(document)

print(f"Signature: {signature}")
print()

if X.verify(signature, document, X.get_public()):
    print("The verification was successful.")
else:
    print("The verification was not successful.")
