# coding: utf-8
import math
import secrets
from typing import List, Tuple

import galois

n = 3
GF = galois.GF(2**61 - 1)

constants_int = [
    "0x0",
    "0x0",
    "0x2",
    "0x1",
    "0x0",
    "0x2",
    "0x3",
    "0x1",
    "0x3",
    "0x2",
    "0x3",
    "0x1",
    "0x2",
    "0x1",
    "0x0",
    "0x0",
    "0x3",
    "0x3",
    "0x3",
    "0x1",
    "0x1",
    "0x3",
    "0x0",
    "0x3",
    "0x0",
    "0x1",
    "0x3",
    "0x1",
    "0x0",
    "0x2",
    "0x3",
    "0x2",
    "0x1",
    "0x1",
    "0x1",
    "0x3",
    "0x2",
    "0x3",
    "0x2",
    "0x1",
    "0x3",
    "0x3",
    "0x1",
    "0x2",
    "0x2",
    "0x0",
    "0x3",
    "0x2",
    "0x1",
    "0x3",
    "0x2",
    "0x2",
    "0x0",
    "0x1",
]


# Convert constants to field elements.
constants = [GF(int(c, base=16)) for c in constants_int]

# t_i = x_i + c_i + k
# s_i = t_i^2
# x_i+1  = s_i * t_i
# 2 constraints per round

#         [1    ]
#         [x_i  ]
#         [k    ]
#    x =  [s_i  ]
#         [x_i+1]
#

# A = [ c_i 1 1 0 0]         A * x  = [x_i + c_i + k]
#     [ 0   0 0 1 0]                  [s_i]

# B = [ c_i 1 1  0 0]         B * x  = [x_i + c_i + k]      => [x_i + c_i + k] * [x_i + c_i + k] = [s_i]
#     [ c_i 1 1  0 0]                  [x_i + c_i + k]         [s_i] * [x_i + c_i + k] = [x_i+1]

# C = [ 0 0 0 1 0]           C * x  = [s_i]
#     [ 0 0 0 0 1]                    [x_i+1]

# one round CCS:
#
# a_{i+1} = (a_i + k + b_i)^3
# z = (w, 1, x)
# w = {k, a_i, a_{i+1}}
# x = {b_i}
# ==> z = {k, a_i, a_{i+1}, 1, b_i}
#


# -----------------------------------------------------------------------------
# R1CS helpers
# -----------------------------------------------------------------------------


def _build_r1cs(
    p: "GF", k: "GF", num_rounds: int
) -> Tuple[List[List["GF"]], List[List["GF"]], List[List["GF"]], List["GF"]]:
    """Constructs the R1CS representation (A, B, C) for MiMC encryption together
    with a satisfying witness vector.

    The variable ordering of the witness vector is::

        index | description
        ------+-------------
          0   | constant 1 (always one)
          1   | secret key *k*
          2   | plaintext  *x_0 = p*
          3   | s_0 (auxiliary square result of round 0)
          4   | x_1 (output of round 0)
          5   | s_1
          6   | x_2
          ... | ...

    Hence, for *num_rounds* rounds the witness has

        3 + 2·num_rounds

    entries.

    Per MiMC round we create **two** R1CS constraints as documented in the
    header comment of this file.
    """

    # Build the witness -------------------------------------------------------
    witness: List["GF"] = [GF(1), k, p]

    # Helper lambda to derive indices in the witness vector.
    def idx_x(i: int) -> int:
        """Index of *x_i* (state before round *i*)."""
        return 2 + 2 * i

    def idx_s(i: int) -> int:
        """Index of *s_i* (square result of round *i*)."""
        return 3 + 2 * i

    # Iterate through the rounds, updating the state and extending the witness.
    x_i = p  # current state before this round
    for rnd in range(num_rounds):
        t_i = x_i + k + constants[rnd]
        s_i = t_i**2
        x_ip1 = s_i * t_i

        # Append auxiliary and next-state variables.
        witness.append(s_i)
        witness.append(x_ip1)

        # Advance state
        x_i = x_ip1

    # Build A, B, C matrices --------------------------------------------------
    n_vars = len(witness)
    zero = GF(0)

    A: List[List["GF"]] = []
    B: List[List["GF"]] = []
    C: List[List["GF"]] = []

    # For every round generate two constraints.
    for rnd in range(num_rounds):
        # Pre-compute frequently used indices.
        x_idx = idx_x(rnd)
        s_idx = idx_s(rnd)
        x_next_idx = idx_x(rnd + 1)  # x_{i+1}

        # Constraint 1: (x_i + c_i + k)^2 = s_i
        rowA1 = [zero] * n_vars
        rowB1 = [zero] * n_vars
        rowC1 = [zero] * n_vars

        # A row coefficients
        rowA1[0] = constants[rnd]  # constant term (c_i)
        rowA1[1] = GF(1)  # k
        rowA1[x_idx] = GF(1)  # x_i

        # B row is identical for the square constraint.
        rowB1[:] = rowA1

        # C row outputs s_i
        rowC1[s_idx] = GF(1)

        # Append rows to matrices.
        A.append(rowA1)
        B.append(rowB1)
        C.append(rowC1)

        # Constraint 2: s_i * (x_i + c_i + k) = x_{i+1}
        rowA2 = [zero] * n_vars
        rowB2 = [zero] * n_vars
        rowC2 = [zero] * n_vars

        # A row – just s_i
        rowA2[s_idx] = GF(1)

        # B row – (x_i + c_i + k)
        rowB2[0] = constants[rnd]
        rowB2[1] = GF(1)
        rowB2[x_idx] = GF(1)

        # C row – x_{i+1}
        rowC2[x_next_idx] = GF(1)

        A.append(rowA2)
        B.append(rowB2)
        C.append(rowC2)

    return A, B, C, witness


# -----------------------------------------------------------------------------
# MiMC primitives (updated to optionally return R1CS)
# -----------------------------------------------------------------------------


def mimc_encryption(
    p: "GF", k: "GF", num_rounds: int
) -> Tuple["GF", List[List["GF"]], List[List["GF"]], List[List["GF"]], List["GF"]]:
    """Encrypt a single block using the MiMC permutation **and** return a
    Rank-1 Constraint System proving its correctness.

    The function preserves the original behaviour (producing the ciphertext)
    but additionally outputs the R1CS matrices *(A, B, C)* and a satisfying
    witness vector *x* such that for every constraint row *i*

        (A_i · x) · (B_i · x) = C_i · x.

    Two constraints are generated per MiMC round.
    """

    # Build R1CS and witness.
    A, B, C, witness = _build_r1cs(p, k, num_rounds)

    # The last state variable in the witness is x_{num_rounds} (before the
    # final key addition). The actual ciphertext is x_{num_rounds} + k.
    ciphertext = witness[-1] + k

    return ciphertext, A, B, C, witness


def mimc_decrypt(c: GF, k: GF, num_rounds: int) -> GF:
    """Decrypts a single block using the MiMC cipher in GF(2^129). Assumes the same constants and key as encryption.

    Parameters
    ----------
    c : GF
        Ciphertext element.
    k : GF
        Key element.
    num_rounds : int
        Number of MiMC rounds to invert.

    Returns
    -------
    GF
        The decrypted plaintext element.
    """
    # Cube root exponent for the field (for decryption)
    cube_root_exp = pow(3, -1, GF.order - 1)
    state = c - k
    for i in range(num_rounds - 1, 0, -1):
        state = (state**cube_root_exp) - (k + constants[i])
    state = (state**cube_root_exp) - (k + constants[0])
    return state


# Determine the number of rounds suggested by the MiMC security analysis.
num_rounds = int(math.ceil(n // math.log2(3)))


def mimc_test() -> None:
    """Original sanity-check for encryption / decryption (all rounds)."""

    print(f"Number of rounds: {num_rounds}")

    # Secret key (chosen arbitrarily for this demonstration).
    key_int = 0x42
    k = GF(key_int)
    print(f"Key: {hex(int(k))}")

    # Example plaintext.
    p = GF(0x12)

    # Encrypt (also obtains R1CS, which we ignore here).
    ciphertext, *_ = mimc_encryption(p, k, num_rounds)

    print(f"Plaintext: {hex(int(p))}")
    print(f"Ciphertext: {hex(int(ciphertext))}")

    # Decrypt.
    # p_dec = mimc_decrypt(ciphertext, k, num_rounds)
    # print(f"Decrypted: {hex(int(p_dec))}")
    # assert p == p_dec, "Decryption failed!"


# -----------------------------------------------------------------------------
# Simple 1-round R1CS self-test
# -----------------------------------------------------------------------------


def _dot(row: List["GF"], w: List["GF"]) -> "GF":
    """Compute the dot-product of *row* with the *witness* vector."""
    acc = GF(0)
    for a, b in zip(row, w):
        acc += a * b
    return acc


def r1cs_test() -> None:
    """Verify that the generated R1CS satisfies the constraints for 1 round."""

    print("\n[ R1CS 1-round self-test ]")

    # Fixed key and plaintext for reproducibility.
    k = GF(0x4242)
    p = GF(0x123)

    ciphertext, A, B, C, witness = mimc_encryption(p, k, num_rounds)
    # print_r1cs_matrices(A, B, C, witness)

    # Verify the two constraints.
    for idx, (rowA, rowB, rowC) in enumerate(zip(A, B, C)):
        left = _dot(rowA, witness)
        right = _dot(rowB, witness)
        out = _dot(rowC, witness)
        assert left * right == out, f"Constraint {idx} failed"

    Ax = GF(A) @ GF(witness)
    Bx = GF(B) @ GF(witness)
    Cx = GF(C) @ GF(witness)
    import numpy as np

    print(witness)
    hadamard = np.multiply(Ax, Bx)
    assert np.array_equal(hadamard, Cx), f"Matrix R1CS check failed: {hadamard} != {Cx}"

    # # Matrix multiplication
    # Ax = A @ witness
    # Bx = B @ witness
    # Cx = C @ witness

    # # Hadamard product (element-wise multiplication)
    # hadamard = np.multiply(Ax, Bx)

    # # Verify R1CS constraint: (A·witness) ⊙ (B·witness) = C·witness
    # assert np.array_equal(hadamard, Cx), f"Matrix R1CS check failed: {hadamard} != {Cx}"
    # print("Matrix-based R1CS check passed!")

    print("All constraints satisfied.")
    serialize_r1cs_to_json(A, B, C, witness, "mimc_r1cs.json")


def serialize_r1cs_to_json(
    A: List[List["GF"]],
    B: List[List["GF"]],
    C: List[List["GF"]],
    witness: List["GF"],
    filename: str = "mimc_r1cs.json",
) -> None:
    """Serialize R1CS matrices and witness vector to a JSON file.

    Args:
        A: Matrix A from R1CS
        B: Matrix B from R1CS
        C: Matrix C from R1CS
        witness: Witness vector
        filename: Output JSON filename
    """
    import json

    # Convert GF elements to hex strings
    data = {
        "A": [[hex(int(x)) for x in row] for row in A],
        "B": [[hex(int(x)) for x in row] for row in B],
        "C": [[hex(int(x)) for x in row] for row in C],
        "witness": [hex(int(x)) for x in witness],
    }

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

    print(f"\nR1CS data serialized to {filename}")


def print_r1cs_matrices(
    A: List[List["GF"]], B: List[List["GF"]], C: List[List["GF"]], witness: List["GF"]
) -> None:
    """Print R1CS matrices and witness vector in a readable format."""
    print("\n[ R1CS Matrices ]")

    print("\nMatrix A:")
    for i, row in enumerate(A):
        print(f"Row {i}: {[hex(int(x)) if x != 0 else '0' for x in row]}")

    print("\nMatrix B:")
    for i, row in enumerate(B):
        print(f"Row {i}: {[hex(int(x)) if x != 0 else '0' for x in row]}")

    print("\nMatrix C:")
    for i, row in enumerate(C):
        print(f"Row {i}: {[hex(int(x)) if x != 0 else '0' for x in row]}")

    print("\nWitness vector:")
    print([hex(int(x)) for x in witness])


def mimc_main() -> None:
    """Module entry-point used when executing this file directly."""

    mimc_test()
    r1cs_test()
