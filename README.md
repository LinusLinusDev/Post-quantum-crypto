# Post Quantum Cryptography

---

### lattice_based.py

---

Implementation of the GHH-Cryptosystem using the HNF as public key:

Mapping an n-dimensional vector of integers onto a given lattice and encoding this point by adding an error.

Recover the point and thus the vector using Babai's rounding procedure.

Calculate the worst possible basis of the given lattice using the Nemhauser/Wolsey algorithm to obtain the Hermite normal form in its lower triangular form.

---

### multivariat.py

---

Implementation of the Unbalanced Oil and Vinegar (UOV) Signature Sheme
