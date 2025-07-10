# coding: utf-8
import math
import secrets
from typing import List, Tuple

import galois

n = 129
GF = galois.GF(4787605948707450321761805915146316350821882368518086721537)

constants_int = [
    0x000000000000000000000000000000000,
    0x06382246D63A75DB1F522FA4E05CC0657,
    0x1081B8FFC1F376F1D8B798232E330ED56,
    0x0CE5D22BD8F8D1D3B262296B92F588787,
    0x1A426AA322E6B4DA51F9DBA2BB404DF1C,
    0x03C6946A088909908BC4471593028E7FB,
    0x11C2EDF7C3929DA2A3A9F3B3B1FC8EB53,
    0x10B0157C0FA825D7E9965ACE273294881,
    0x0AF67AC1897FF10C25501933D6C930E94,
    0x00BBBC3C4BBD5A54A0E66D9F0F8A8BE79,
    0x0CDD6A5BE343BEC16D6183F7BADD6D5DB,
    0x08B843AF944F96E199C9DAC93A8AF2888,
    0x18EA7BD63B43D8C85512B7669D9927518,
    0x16D13F69CFF104E013754C2FE96F5CB92,
    0x168AA09F0672C5E465461AECDD6684B6A,
    0x1FE480707DCA72BFCCB1A77CD4EDB67A6,
    0x015CEE01DACF5C9E66FE6BDFA93189454,
    0x16C891BC68270F24998606CFBDCAA5611,
    0x0CF4FD57AD35BF06D242A9FCB4972075B,
    0x1170A49DB9A62896C1B59840E33FF427D,
    0x0D675BD0083F069681057099A88261932,
    0x1E35E81A5F3B88568CD18456936BE6721,
    0x078A0AC6891B377B855F7E715144133B8,
    0x1664154B3DE279C0D7FC6DD5C1EC4ADA8,
    0x01860E08D1D0CCEB6E0D068D753AFB5B8,
    0x0DDBC94E9E569546044A616AB8F462FE8,
    0x06A5AE6C7C22904E4BAF44848294F82C9,
    0x05ADF151652A16F0E3546A701953DE05B,
    0x1590219323779A7756933E1D43B092865,
    0x0F09EA94A42CF3D8182EC81FE19B4A16B,
    0x14A00ADFAA82483DB455CED5B42588C23,
    0x0F9B2C70CABEA07F738F6C4C2FDF05271,
    0x177D81B15626B3D9705A84BA498BED335,
    0x010E37AD1A843D68C93BD4D12A5FF777E,
    0x1EC5252E4CE05E567C0ED58E392F82DF9,
    0x0C51762CDB2A41B86FADCA23ADE46EC3C,
    0x1BD09CBC69B6AA86C79D4B56E06C65DBF,
    0x0472CF4DF04D10A8764DCE39C758AC89F,
    0x0FA41F59A323DFFDC05BDCA384BA65007,
    0x1056526AA50107101EB34D698DBB1507B,
    0x10F3E7E57D9451BB1C36A7DB077D623B2,
    0x0C843A8A873AA33444962D64B243CA1F1,
    0x1A06A6D9E59A6F17F829845829FCC51A5,
    0x078D4F7FA105E1F396B2BDEC55D07E96A,
    0x1336A5EB8A15D2AE237E6605B4A4D5E7E,
    0x084BD7DDFC3F58851E405EEE24B31E0DA,
    0x196F9C6845D9ABC8B17815E4EFE43EA61,
    0x17E060968262D38D5B12BE87B0DDCA0FE,
    0x13E5CA95C7826B284615893FB6B6615F2,
    0x11AAFE1DEF56E71FCF9FFD1F535472262,
    0x0538E78611F47E30797CBFA5EB217D9FD,
    0x1ADB78502384A7A093D4CB5423EB98DFE,
    0x03F564C552C72F1F615660913F31AE19E,
    0x05383AEC281ADD2A5E61FC16E0915E9AB,
    0x17B7D2F156B797BA3BCD04F74970A3698,
    0x16BCB475655EECE3A2F8DECD844F65550,
    0x0DF93054F75B723AC4E2CE48D00CF37DD,
    0x190B65B81EF953C92AD0B5A15C533824C,
    0x026E03D771818ACFAD02DD38D3D5AD6D7,
    0x088946EE4840404BF1FE6EF874751680F,
    0x1A2164A4A31C13D1A0FE4D86B8A5A8F0C,
    0x0C7B325E4AECB36F489A24A31277C18AC,
    0x1A1B145F688B87D5E5926BD19D70858F5,
    0x166976D9031782C3A733897C19EADF660,
    0x0CB0DE4A36207611A580A97D94A99708E,
    0x0453B6E0F8FB6B59D38B466B9B4210B4E,
    0x13BBDD7CDA3D39A2BC6391929699B1D02,
    0x05449506BB6FA430FF999D13EF9187631,
    0x152C672A79FB3A4DE06BF22E9F8A6F7DD,
    0x09CF98E8DB80E7EC38C662CF0BD84DD49,
    0x07C42B47224719B2E7D6416E7AECE843D,
    0x18D14C8C96531EE939835090C92A79A08,
    0x189AC9A8952DAFA06B3FAD1ABE9CF37A8,
    0x0B382A9F685108884C841CFCDD4E7C065,
    0x0263C639FAE4BEE461BC66BE8FED407F7,
    0x118BBB5A626F4130A3246BF144DDEBA6F,
    0x11C7F739620FC72BA7112461FB96BCEF2,
    0x1CE202833557D1E76AF8A03CF4E1FCCF7,
    0x0A474673A25C26E1C18AEAB2015ADDA20,
    0x15C9722C814B888297FCC8C2A096A8730,
    0x1DE01E75FA74625E8F0D8231A510C88DC,
    0x12057179D8D7584FCDFB1C7302988550A,
]

# Convert constants to field elements.
constants = [GF(c) for c in constants_int]

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


# Cube root exponent for the field (for decryption)
cube_root_exp = pow(3, -1, GF.order - 1)


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
    state = c - k
    for i in range(num_rounds - 1, 0, -1):
        state = (state**cube_root_exp) - (k + constants[i])
    state = (state**cube_root_exp) - (k + constants[0])
    return state


# Determine the number of rounds suggested by the MiMC security analysis.
num_rounds = int(math.ceil(n / math.log2(3)))


def mimc_test() -> None:
    """Original sanity-check for encryption / decryption (all rounds)."""

    print(f"Number of rounds: {num_rounds}")

    # Secret key (chosen arbitrarily for this demonstration).
    key_int = 0x42424242424242424242424242424242
    k = GF(key_int)
    print(f"Key: {hex(int(k))}")

    # Example plaintext.
    p = GF(0x123123123123)

    # Encrypt (also obtains R1CS, which we ignore here).
    ciphertext, *_ = mimc_encryption(p, k, num_rounds)

    print(f"Plaintext: {hex(int(p))}")
    print(f"Ciphertext: {hex(int(ciphertext))}")

    # Decrypt.
    p_dec = mimc_decrypt(ciphertext, k, num_rounds)
    print(f"Decrypted: {hex(int(p_dec))}")
    assert p == p_dec, "Decryption failed!"


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
    k = GF(0x42424242424242424242424242424242)
    p = GF(0x123123123123)

    ciphertext, A, B, C, witness = mimc_encryption(p, k, num_rounds)
    # print_r1cs_matrices(A, B, C, witness)

    # Verify the two constraints.
    for idx, (rowA, rowB, rowC) in enumerate(zip(A, B, C)):
        left = _dot(rowA, witness)
        right = _dot(rowB, witness)
        out = _dot(rowC, witness)
        assert left * right == out, f"Constraint {idx} failed"

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
