# The MiMc R1CS constraint and witness is calculated in the Mersenne 61 prime field F_q, where q = 2^61 - 1.
# They are stored in json format, with the following structure:
# {
#     "A": [[int, int, ...], [int, int, ...], ...], # constraint matrix A
#     "B": [[int, int, ...], [int, int, ...], ...], # constraint matrix B
#     "C": [[int, int, ...], [int, int, ...], ...], # constraint matrix C
#     "witness": [int, int, ...] # witness vector
# }
# transform_r1cs_to_low_norm is a generic function that can be used to transform the R1CS constraint and witness to low norm.
# It takes two arguments:
# - q: the field size of the R1CS constraint and witness.
# - norm_bound: the norm bound of the R1CS constraint and witness. constraints and witness are bounded by norm_bound.
# - r1cs_json_path: the path to the json file that contains the R1CS constraint and witness.
# It returns the transformed R1CS constraint and witness.
# The test case is to bound the norm of the R1CS constraint and witness to 2^54.
# This means that all values in the R1CS constraint and the witness must be less than 2^54,
# In other words, given A,B,C, w, the algorithm should find a new A', B', C', w' such that:
# 1. A' * w' \hardamand B' * w' =  C' * w'
# 2. The norm of A', B', C' is less than norm_bound.

import math
import json
import numpy as np


# -----------------------------------------------------------------------------
# A minimal builder to accumulate the low-norm R1CS we are constructing.
# It purposefully uses a *sparse* representation (dict per row) so we do not
# blow up memory on big circuits.
# -----------------------------------------------------------------------------


class R1CSBuilder:
    """Light-weight R1CS builder.

    The constant 1 is assigned index 0 in the witness vector so that we can add
    constant terms the usual R1CS way (coefficient on variable 0).
    """

    def __init__(self, q: int, d: int):
        self.q = q  # field modulus
        self.d = d  # limb base (= norm_bound)

        #   witness[0] == 1  is the conventional constant wire.
        self.witness = [1]

        # A′, B′, C′ rows as sparsely-filled dicts {var_idx → coeff}
        self.A_rows: list[dict[int, int]] = []
        self.B_rows: list[dict[int, int]] = []
        self.C_rows: list[dict[int, int]] = []

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

    def new_var(self, value: int) -> int:
        """Append *value* to the witness vector and return its index."""

        assert 0 <= value < self.q, "value not in field"
        idx = len(self.witness)
        self.witness.append(int(value))
        return idx

    # Alias – we sometimes want to emphasise constants
    def new_const(self, value: int) -> int:  # noqa: D401
        return self.new_var(value)

    @property
    def const_one_idx(self) -> int:  # noqa: D401
        return 0

    def get_val(self, idx: int) -> int:
        return self.witness[idx]

    # ------------------------------------------------------------------
    # Constraint emission
    # ------------------------------------------------------------------

    def add_constraint(
        self,
        a_row: dict[int, int],
        b_row: dict[int, int],
        c_row: dict[int, int],
    ) -> None:
        """Append a row A·w  ∘  B·w  =  C·w ."""

        self.A_rows.append(a_row)
        self.B_rows.append(b_row)
        self.C_rows.append(c_row)

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------

    def to_matrices(self):
        """Return (A',B',C') as list-of-dict rows and the witness list."""

        return self.A_rows, self.B_rows, self.C_rows, self.witness


# -----------------------------------------------------------------------------
# Gadgets (all work with *variable-length* limb arrays)
# -----------------------------------------------------------------------------


def add_range_check_constraints(builder: R1CSBuilder, limb_indices: list[int]):
    """Very cheap range check: host-side assert + identity row."""

    for idx in limb_indices:
        val = builder.get_val(idx)
        assert val < builder.d, "limb exceeds bound"

        # Identity constraint: limb * 1 = limb  ⇢  {idx:1}·w ∘ {0:1}·w = {idx:1}·w
        builder.add_constraint({idx: 1}, {builder.const_one_idx: 1}, {idx: 1})


def add_multilimb_addition(
    builder: R1CSBuilder, a_limbs: list[int], b_limbs: list[int]
) -> list[int]:
    """Return indices of the limbs of (A + B) in base-d, variable-length."""

    # Pad the shorter list with *constant 0* limbs
    zero_idx = builder.new_const(0)
    max_len = max(len(a_limbs), len(b_limbs))
    a = a_limbs + [zero_idx] * (max_len - len(a_limbs))
    b = b_limbs + [zero_idx] * (max_len - len(b_limbs))

    carry_idx = zero_idx
    result: list[int] = []

    for i in range(max_len):
        a_i, b_i = a[i], b[i]

        # ----------------------------------------------------------
        # (1) raw limb sum s_i = a_i + b_i + carry_{i-1}
        # ----------------------------------------------------------
        s_val = (
            builder.get_val(a_i) + builder.get_val(b_i) + builder.get_val(carry_idx)
        ) % builder.q
        s_idx = builder.new_var(s_val)
        builder.add_constraint(
            {a_i: 1, b_i: 1, carry_idx: 1, s_idx: -1},
            {builder.const_one_idx: 1},
            {},
        )

        # ----------------------------------------------------------
        # (2) Decompose s_i into (p_i  +  carry_i * d)
        # ----------------------------------------------------------
        p_val = s_val % builder.d
        c_val = s_val // builder.d
        p_idx = builder.new_var(p_val)
        carry_idx = builder.new_var(c_val)

        # constraint  c*d + p - s = 0  (in A side)
        builder.add_constraint(
            {carry_idx: builder.d, p_idx: 1, s_idx: -1},
            {builder.const_one_idx: 1},
            {},
        )

        # range-check the two fresh limbs
        add_range_check_constraints(builder, [p_idx, carry_idx])

        result.append(p_idx)

    # Final carry (if non-zero)
    if builder.get_val(carry_idx):
        result.append(carry_idx)

    return result


def _recompose(builder: R1CSBuilder, limb_indices: list[int]) -> int:
    """Utility: turn base-d limbs back into a single integer."""

    acc = 0
    for i, idx in enumerate(limb_indices):
        acc += builder.get_val(idx) * (builder.d**i)
    return acc % builder.q


def add_multilimb_multiplication(
    builder: R1CSBuilder, a_limbs: list[int], b_limbs: list[int]
) -> list[int]:
    """Multiply two multi-limb numbers and return limb indices of the product."""

    # Strategy for prototype simplicity: 1) compute the integer product in
    # Python, 2) decompose, 3) introduce all needed constraints linking the
    # operands and the result via a single R1CS multiplication + equality.

    a_val = _recompose(builder, a_limbs)
    b_val = _recompose(builder, b_limbs)
    prod_val = (a_val * b_val) % builder.q

    prod_limbs_vals = decompose(prod_val, builder.d)
    prod_idx = [builder.new_var(v) for v in prod_limbs_vals]

    # Range-check result limbs.
    add_range_check_constraints(builder, prod_idx)

    # R1CS enforcement: introduce *single* variable holding reconstructed prod
    prod_full_idx = builder.new_var(prod_val)
    # prod_full == a * b   via constraint   a_full · b_full = prod_full
    a_full_idx = builder.new_var(a_val)
    b_full_idx = builder.new_var(b_val)
    builder.add_constraint({a_full_idx: 1}, {b_full_idx: 1}, {prod_full_idx: 1})

    # Link prod_full to limbs:  Σ limb_i * d^i  == prod_full.
    lhs = {prod_full_idx: -1}
    for i, limb in enumerate(prod_idx):
        coeff = builder.d**i
        lhs[limb] = coeff
    builder.add_constraint(lhs, {builder.const_one_idx: 1}, {})

    return prod_idx


def add_multilimb_dot_product(
    builder: R1CSBuilder,
    coeff_dict: dict[int, list[int]],
    witness_map: dict[int, list[int]],
) -> list[int]:
    """Return limbs of Σ_j coeff_j * w_j ."""

    # Accumulate partial products as limbs via the addition gadget.

    result_limbs: list[int] | None = None
    for var_idx, coeff_limbs_vals in coeff_dict.items():
        # (1) make constant variables for coefficient limbs
        coeff_idx = [builder.new_const(v) for v in coeff_limbs_vals]

        # (2) fetch witness limbs
        w_limbs = witness_map[var_idx]

        # (3) multiply ⇒ term limbs
        term_limbs = add_multilimb_multiplication(builder, coeff_idx, w_limbs)

        # (4) accumulate
        if result_limbs is None:
            result_limbs = term_limbs
        else:
            result_limbs = add_multilimb_addition(builder, result_limbs, term_limbs)

    # Edge case: if no coeff ≥ d, the dot product is 0
    if result_limbs is None:
        zero_idx = builder.new_const(0)
        result_limbs = [zero_idx]

    return result_limbs


def add_multilimb_multiplication_constraint(
    builder: R1CSBuilder,
    L_limbs: list[int],
    R_limbs: list[int],
    O_limbs: list[int],
):
    """Enforce that the multi-limb product of L and R equals O."""

    prod_limbs = add_multilimb_multiplication(builder, L_limbs, R_limbs)

    # subtract O limbs: prod − O == 0  limb-wise
    max_len = max(len(prod_limbs), len(O_limbs))
    zero_idx = builder.new_const(0)

    prod = prod_limbs + [zero_idx] * (max_len - len(prod_limbs))
    O = O_limbs + [zero_idx] * (max_len - len(O_limbs))

    diff_limbs = add_multilimb_addition(
        builder,
        prod,
        [
            builder.new_const((-builder.get_val(idx)) % builder.q) for idx in O
        ],  # negate O
    )

    # Enforce every limb of diff is 0
    for idx in diff_limbs:
        builder.add_constraint({idx: 1}, {builder.const_one_idx: 1}, {})


def transform_r1cs_to_low_norm(A, B, C, w, q, norm_bound):
    """Convert (A,B,C,w) into an equivalent low-norm R1CS."""

    # ------------------------------------------------------------------
    # 1. Builder & bookkeeping
    # ------------------------------------------------------------------
    builder = R1CSBuilder(q, norm_bound)
    witness_map: dict[int, list[int]] = {}

    # ------------------------------------------------------------------
    # 2. Witness transformation – variable-length limb encoding
    # ------------------------------------------------------------------
    for j, w_j in enumerate(w):
        if w_j >= norm_bound:
            limb_vals = decompose(w_j, norm_bound)
        else:
            limb_vals = [w_j]

        # Sanity checks
        assert w_j == sum(l * norm_bound**i for i, l in enumerate(limb_vals))
        assert all(l < norm_bound for l in limb_vals)

        # Materialise variables in the builder & range-check them
        limb_idx = [builder.new_var(v) for v in limb_vals]
        add_range_check_constraints(builder, limb_idx)

        witness_map[j] = limb_idx

    # ------------------------------------------------------------------
    # 3. Constraint transformation – row by row
    # ------------------------------------------------------------------
    num_constraints = A.shape[0]
    assert num_constraints == B.shape[0] == C.shape[0]

    for k in range(num_constraints):
        # Decompose *coefficients* ≥ bound (they are constants, so we keep them
        # as values for now).  Coefficients < bound stay single limb.
        a_k_limbs = {
            j: ([A[k, j]] if A[k, j] < norm_bound else decompose(A[k, j], norm_bound))
            for j in A[k].nonzero()[0]
        }
        b_k_limbs = {
            j: ([B[k, j]] if B[k, j] < norm_bound else decompose(B[k, j], norm_bound))
            for j in B[k].nonzero()[0]
        }
        c_k_limbs = {
            j: ([C[k, j]] if C[k, j] < norm_bound else decompose(C[k, j], norm_bound))
            for j in C[k].nonzero()[0]
        }

        # Build limb representation of the dot products
        L_limbs = add_multilimb_dot_product(builder, a_k_limbs, witness_map)
        R_limbs = add_multilimb_dot_product(builder, b_k_limbs, witness_map)
        O_limbs = add_multilimb_dot_product(builder, c_k_limbs, witness_map)

        # Enforce multiplication equality
        add_multilimb_multiplication_constraint(builder, L_limbs, R_limbs, O_limbs)

    # ------------------------------------------------------------------
    # 4. Export the new system
    # ------------------------------------------------------------------
    return builder.to_matrices()


def decompose(value, base):
    """Return limbs of ``value`` in the given ``base`` until the value is exhausted.

    At least one limb is returned, so ``0`` maps to ``[0]``.
    """
    if value == 0:
        return [0]

    limbs = []
    temp_val = value
    while temp_val > 0:
        limbs.append(temp_val % base)
        temp_val //= base
    return limbs


# -----------------------------------------------------------------------------
# R1CS verification helper
# -----------------------------------------------------------------------------


def _dot_row(row, witness, q):
    """Compute Σ coeff_i * w_i (mod q) for a *single* R1CS row.

    The helper seamlessly supports three common encodings for an R1CS row:

    1. Sparse dict:          ``{var_idx: coeff, ...}``
    2. NumPy / list dense:   ``np.ndarray`` or ``list`` of coefficients

    The dense encoding *may* be large in practice.  For efficiency we skip
    zero coefficients when the container exposes ``nonzero`` (NumPy).
    """
    # Sparse representation -------------------------------------------------
    if isinstance(row, dict):
        return sum(int(coeff) * int(witness[idx]) for idx, coeff in row.items()) % q

    # NumPy ndarray ---------------------------------------------------------
    import numpy as _np  # local import to avoid mandatory dependency at file import time

    if isinstance(row, _np.ndarray):
        # Quickly extract the non-zero indices to avoid a full scan of large rows.
        nz_idx = row.nonzero()[0]
        return sum(int(row[i]) * int(witness[i]) for i in nz_idx) % q

    # Generic sequence (e.g. Python list) ----------------------------------
    return sum(int(coeff) * int(witness[i]) for i, coeff in enumerate(row) if coeff) % q


def verify_r1cs_constraints(A_rows, B_rows, C_rows, witness, q):
    """Verify that (A·w) ⊙ (B·w) = C·w.

    The function accepts *either* sparse (dict-per-row) or dense (NumPy array
    or Python list) encodings for the three matrices.  The witness vector may
    likewise be any sequence supporting integer indexing.
    """

    for k, (a_row, b_row, c_row) in enumerate(zip(A_rows, B_rows, C_rows)):
        left = _dot_row(a_row, witness, q)
        right = _dot_row(b_row, witness, q)
        out = _dot_row(c_row, witness, q)

        if (left * right) % q != out:
            raise AssertionError(
                f"Constraint {k} failed: {(left * right) % q} != {out}"
            )


if __name__ == "__main__":

    q = 2**61 - 1
    norm_bound = 2**54

    with open("mimc_r1cs.json", "r") as f:
        data = json.load(f)
        A = np.array([[int(x, 16) for x in row] for row in data["A"]])
        B = np.array([[int(x, 16) for x in row] for row in data["B"]])
        C = np.array([[int(x, 16) for x in row] for row in data["C"]])
        witness = np.array([int(x, 16) for x in data["witness"]])

    verify_r1cs_constraints(A, B, C, witness, q)
    A_p, B_p, C_p, w_p = transform_r1cs_to_low_norm(A, B, C, witness, q, norm_bound)
    verify_r1cs_constraints(A_p, B_p, C_p, w_p, q)

    print("All constraints satisfied!")
    print("#constraints:", len(A_p))
    print("#witness variables:", len(w_p))
    # simple multiplication test

    witness = np.array([1, 2**55, 3, 3 * 2**55])
    A = np.array([[0, 1, 0, 0]])
    B = np.array([[0, 0, 1, 0]])
    C = np.array([[0, 0, 0, 1]])
    verify_r1cs_constraints(A, B, C, witness, q)
    A_p, B_p, C_p, w_p = transform_r1cs_to_low_norm(A, B, C, witness, q, norm_bound)
    verify_r1cs_constraints(
        A_p, B_p, C_p, w_p, q
    )  # verify the transformed R1CS constraints are satisfied

    print("All constraints satisfied!")
    print("#constraints:", len(A_p))
    print("#witness variables:", len(w_p))

    # ------------------------------------------------------------------
    # Additional test 1 – large witness values
    # ------------------------------------------------------------------

    a = 2**60 - 123  # ≥ norm_bound, triggers decomposition
    b = 2**60 - 45678  # ≥ norm_bound, triggers decomposition
    prod = (a * b) % q

    witness2 = np.array([1, a, b, prod], dtype=object)
    A2 = np.array([[0, 1, 0, 0]], dtype=object)
    B2 = np.array([[0, 0, 1, 0]], dtype=object)
    C2 = np.array([[0, 0, 0, 1]], dtype=object)

    verify_r1cs_constraints(A2, B2, C2, witness2, q)
    A2_p, B2_p, C2_p, w2_p = transform_r1cs_to_low_norm(
        A2, B2, C2, witness2, q, norm_bound
    )
    verify_r1cs_constraints(A2_p, B2_p, C2_p, w2_p, q)

    print("\nLarge witness test passed!")
    print("#constraints:", len(A2_p))
    print("#witness variables:", len(w2_p))

    # ------------------------------------------------------------------
    # Additional test 2 – large coefficients in the matrices
    # ------------------------------------------------------------------
    big_coeff = 2**60 - 321  # ≥ norm_bound, appears as coefficient
    x = 2**60 - 987  # witness value ≥ norm_bound
    y = (big_coeff * x) % q

    witness3 = np.array([1, x, y], dtype=object)

    # A row applies big_coeff to x, B row multiplies by constant wire 1, C row expects the product
    A3 = np.array([[0, big_coeff, 0]], dtype=object)  # coef on x
    B3 = np.array([[1, 0, 0]], dtype=object)  # constant 1 (wire 0)
    C3 = np.array([[0, 0, 1]], dtype=object)  # coef on y

    verify_r1cs_constraints(A3, B3, C3, witness3, q)
    A3_p, B3_p, C3_p, w3_p = transform_r1cs_to_low_norm(
        A3, B3, C3, witness3, q, norm_bound
    )
    verify_r1cs_constraints(A3_p, B3_p, C3_p, w3_p, q)

    print("\nLarge coefficient test passed!")
    print("#constraints:", len(A3_p))
    print("#witness variables:", len(w3_p))
