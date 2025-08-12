# The MiMc R1CS constraint and witness is calculated in the Mersenne 61 prime field F_q, where q = 2^61 - 1.
# They are stored in json format, with the following structure:
# {
#     "A": [[int, int, ...], [int, int, ...], ...], # constraint matrix A
#     "B": [[int, int, ...], [int, int, ...], ...], # constraint matrix B
#     "C": [[int, int, ...], [int, int, ...], ...], # constraint matrix C
#     "witness": [int, int, ...] # witness vector
# }
# This is a generic function that can be used to transform the R1CS constraint and witness to low norm.
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
# The algorithm is as follows:

# witness = [1, 2**55, 3, 3 * 2**55]
# A = [[0, 1, 0, 0]]
# B = [[0, 0, 1, 0]]
# C = [[0, 0, 0, 1]]

# d=2^54

# 6 * 2**55 = 0*(2^54)^0 + 6*(2^54)^1 + ...
# 0*(2^54)^0 + 6*(2^54)^1 + 7*(2^54)^2...

# 2^55 = [0*d^0, 2*d^1]

# w'=[1, 0, 2, 3, 0, 6]
# witness_map = {0: [0], 1: [1, 2], 2: [3], 3: [4, 5]}

# A' = [[0, 1*d^0, 1*d^1, 0, 0, 0]]
# B' = [[0, 0, 1, 0, 0, 0]]
# C' = [[0, 0, 0, 0, 1*d^0, 1*d^1]]

# but this won't work for arbitrary b since further decomposition would give higher degrees of b
# d^2, d^3,... and this would violate the bound

# But we can always come up with a d' < d, such that a decomposition
# for the largest element in the R1CS is possible without violating the bound d
# i.e. the largest element in the R1CS is decomposed in base d' with
# largest element of this decomposition (d')**k < d


import math
import json


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

        # Pre-populate witness with constant 1 so that index 0 is the
        # usual R1CS constant wire.
        self.witness = [1]

        # A′, B′, C′ rows as sparsely-filled dicts {var_idx → coeff}
        self.A_rows: list[dict[int, int]] = []
        self.B_rows: list[dict[int, int]] = []
        self.C_rows: list[dict[int, int]] = []
        # Shortcut lists so we can refer to “matrix index → row list”.
        self.M = [self.A_rows, self.B_rows, self.C_rows]

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
        i: int,
        row: dict[int, int],
    ) -> None:
        self.M[i].append(row)

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------

    def to_matrices(self):
        """Return (A',B',C') as list-of-dict rows and the witness list."""

        return self.A_rows, self.B_rows, self.C_rows, self.witness


# -----------------------------------------------------------------------------
# Gadgets (all work with *variable-length* limb arrays)
# -----------------------------------------------------------------------------


def add_multilimb_addition(
    builder: R1CSBuilder, a_limbs: list[int], b_limbs: list[int]
) -> list[int]:
    """Create limb variables for the sum of two multi-limb numbers.

    The implementation is *very* light-weight: we recompute the integer
    values off-circuit, decompose the result in base-d and materialise new
    limb variables.  We *do not* yet add carry/overflow correctness
    constraints – this will be done once a proper range-check gadget is in
    place.  To make sure every limb shows up in the final R1CS we emit a
    trivial identity constraint  x · 1 = x  for each of the operands and
    for every freshly created result limb.
    """

    base = builder.d

    # -------- 1. Recompose → integers ------------------------------------
    val_a = _recompose(builder, a_limbs)
    val_b = _recompose(builder, b_limbs)

    # -------- 2. Compute sum and decompose --------------------------------
    sum_val = val_a + val_b
    sum_digits = decompose(sum_val, base)

    # -------- 3. Materialise result limbs ---------------------------------
    sum_limbs_idx = [builder.new_var(v % builder.q) for v in sum_digits]

    # -------- 4. Ensure all limbs participate in at least one constraint --
    for idx in (*a_limbs, *b_limbs, *sum_limbs_idx):
        builder.A_rows.append({idx: 1})
        builder.B_rows.append({builder.const_one_idx: 1})
        builder.C_rows.append({idx: 1})

    return sum_limbs_idx


def add_multilimb_var_multiplication(
    builder: R1CSBuilder,
    x_limbs: list[int],
    y_limbs: list[int],
) -> list[int]:
    """Multiply two multi-limb **variables** and materialise product limbs.

    Very similar philosophy to the constant-by-variable multiplication gadget
    we already have: perform the arithmetic off-circuit, decompose, allocate
    new witness variables for each product digit and add a trivial identity
    constraint to keep the variables alive.  Correctness (carry logic) will be
    added later.
    """

    base = builder.d

    # 1. Recompose to integer values
    val_x = _recompose(builder, x_limbs)
    val_y = _recompose(builder, y_limbs)

    # 2. Compute product and decompose in base-d with one extra limb of
    #    head-room (optional) – the decompose helper already strips leading
    #    zeros, so we just append a 0 limb to guarantee the extra space.
    prod_val = val_x * val_y
    prod_digits = decompose(prod_val, base)
    prod_digits.append(0)  # extra carry limb (value < d)

    # 3. Materialise limbs
    prod_limbs_idx = [builder.new_var(v % builder.q) for v in prod_digits]

    # 4. Identity constraints so that limbs participate

    for idx in (*x_limbs, *y_limbs, *prod_limbs_idx):
        builder.A_rows.append({idx: 1})
        builder.B_rows.append({builder.const_one_idx: 1})
        builder.C_rows.append({idx: 1})

    return prod_limbs_idx


# -----------------------------------------------------------------------------
# Hadamard-product constraint
# -----------------------------------------------------------------------------


def add_multilimb_hadamard_constraint(
    builder: R1CSBuilder,
    A_terms: list[list[int]],
    B_terms: list[list[int]],
    C_terms: list[list[int]],
) -> None:
    """Enforce (A·w) ⊙ (B·w) = C·w on limb arrays.

    * ``A_terms`` / ``B_terms`` / ``C_terms`` are lists where each element is a
      list of limb indices produced by ``add_multilimb_multiplication`` – one
      element per non-zero entry in the original R1CS row.
    """

    # ------------------------------------------------------------------
    # 1. Helper to aggregate many terms via repeated addition
    # ------------------------------------------------------------------
    def _sum_all(terms: list[list[int]]) -> list[int]:
        assert len(terms) > 0, "Must have at least one term"
        cur = terms[0]
        for nxt in terms[1:]:
            cur = add_multilimb_addition(builder, cur, nxt)
        return cur

    A_sum = _sum_all(A_terms)
    B_sum = _sum_all(B_terms)
    C_sum = _sum_all(C_terms)

    # ------------------------------------------------------------------
    # 2. Multiply the summed dot-products
    # ------------------------------------------------------------------
    P_limbs = add_multilimb_var_multiplication(builder, A_sum, B_sum)

    # ------------------------------------------------------------------
    # 3. Pad limb vectors to equal length
    # ------------------------------------------------------------------
    L = max(len(P_limbs), len(C_sum))

    def _pad_with_zeroes(arr: list[int]):
        while len(arr) < L:
            zero_idx = builder.new_var(0)
            # trivial identity so the variable is used
            builder.A_rows.append({zero_idx: 1})
            builder.B_rows.append({builder.const_one_idx: 1})
            builder.C_rows.append({zero_idx: 1})
            arr.append(zero_idx)

    _pad_with_zeroes(P_limbs)
    _pad_with_zeroes(C_sum)

    # ------------------------------------------------------------------
    # 4. Enforce limb-wise equality: P_i · 1 = C_i
    # ------------------------------------------------------------------
    for p_idx, c_idx in zip(P_limbs, C_sum):
        builder.A_rows.append({p_idx: 1})
        builder.B_rows.append({builder.const_one_idx: 1})
        builder.C_rows.append({c_idx: 1})


def _recompose(builder: R1CSBuilder, limb_indices: list[int]) -> int:
    """Utility: turn base-d limbs back into a single integer."""

    acc = 0
    for i, idx in enumerate(limb_indices):
        acc += builder.get_val(idx) * (builder.d**i)
    return acc % builder.q


def add_multilimb_multiplication(
    builder: R1CSBuilder,
    m_idx: int,
    c_limbs: dict[int, list[int]],
    w_limbs_idx: list[int],
) -> list[int]:
    """
    Multiply two multi-limb numbers and add constraints to the builder.
    m_idx is the index of the constraint matrix to add the constraints to.
    c_limbs is a dictionary of variable index to a list of limb values in a row of constraint matrix.
    w_limbs_idx is a list of witness variable indices. these indices can be used to fetch the witness values using builder.get_val.
    The function will add constraints to the builder to enforce the multiplication.
    The function will return a list of limb indices of the product.
    """

    # We only support a *single* term (one wire) at a time.
    assert len(c_limbs) == 1, "c_limbs must contain exactly one entry"

    (var_idx, c_digit_vals) = next(iter(c_limbs.items()))

    # ------------------------------------------------------------------
    # 1. Gather constant-side limbs and witness-side limb values
    # ------------------------------------------------------------------
    base = builder.d
    c_digits = [int(v) for v in c_digit_vals]
    assert all(0 <= v < base for v in c_digits), "Coefficient limbs must be < d"

    w_digits = [builder.get_val(idx) for idx in w_limbs_idx]
    assert all(0 <= v < base for v in w_digits), "Witness limbs must be < d"

    # ------------------------------------------------------------------
    # 2. Compute the product limbs (school-book, to obtain prod_limbs_idx)
    # ------------------------------------------------------------------
    res_len = len(c_digits) + len(w_digits) + 1  # +1 for a possible final carry
    accum: list[int] = [0] * res_len

    for i, c in enumerate(c_digits):
        for j, w in enumerate(w_digits):
            accum[i + j] += c * w  # value < d^2, but we allow big integer here

    # Carry-propagate so every limb < d
    for k in range(len(accum) - 1):
        carry = accum[k] // base
        accum[k] %= base
        accum[k + 1] += carry

    # Strip leading zeros (but keep at least 1 limb)
    while len(accum) > 1 and accum[-1] == 0:
        accum.pop()

    assert all(0 <= v < base for v in accum), "Normalisation failed – limb ≥ d"

    # ------------------------------------------------------------------
    # 3. Materialise the result limbs in the builder (digits of the product)
    # ------------------------------------------------------------------
    prod_limbs_idx: list[int] = [builder.new_var(v % builder.q) for v in accum]

    # ------------------------------------------------------------------
    # 4. NEW: Add correctness constraints tying limbs together using powers of d
    # ------------------------------------------------------------------
    # 4.a  Allocate single-variable representatives of the multi-limb numbers
    #      W_val  – value of the witness-side multiplicand
    #      P_val  – value of the product (mod q)
    #      C_val  – value of the constant-side multiplicand (no variable needed)

    C_val = sum(c * pow(base, i, builder.q) for i, c in enumerate(c_digits)) % builder.q
    W_val_int = sum(w * (base**j) for j, w in enumerate(w_digits))
    W_val = W_val_int % builder.q

    w_val_idx = builder.new_var(W_val)

    P_val = (C_val * W_val) % builder.q
    prod_val_idx = builder.new_var(P_val)

    # 4.b  Constraint (1):  Σ d^j · w_j  =  W_val
    rowA1 = {idx: pow(base, j, builder.q) for j, idx in enumerate(w_limbs_idx)}
    rowB1 = {builder.const_one_idx: 1}
    rowC1 = {w_val_idx: 1}
    builder.A_rows.append(rowA1)
    builder.B_rows.append(rowB1)
    builder.C_rows.append(rowC1)

    # 4.c  Constraint (2):  W_val · C_val  =  P_val
    rowA2 = {w_val_idx: 1}
    rowB2 = {builder.const_one_idx: C_val}
    rowC2 = {prod_val_idx: 1}
    builder.A_rows.append(rowA2)
    builder.B_rows.append(rowB2)
    builder.C_rows.append(rowC2)

    # 4.d  Constraint (3):  Σ d^k · p_k  =  P_val
    rowA3 = {idx: pow(base, k, builder.q) for k, idx in enumerate(prod_limbs_idx)}
    rowB3 = {builder.const_one_idx: 1}
    rowC3 = {prod_val_idx: 1}
    builder.A_rows.append(rowA3)
    builder.B_rows.append(rowB3)
    builder.C_rows.append(rowC3)

    return prod_limbs_idx


def add_multilimb_dot_product(
    builder: R1CSBuilder,
    m_idx: int,
    coeff_dict: dict[int, list[int]],
    witness_map: dict[int, list[int]],
) -> list[list[int]]:
    """Return limbs of Σ_j coeff_j * w_j ."""

    # Accumulate partial products as limbs via the addition gadget.

    result_limbs: list[list[int]] = []
    for var_idx, coeff_limbs_vals in coeff_dict.items():
        # (1) fetch constraint matrix limbs
        c_limbs = {var_idx: coeff_limbs_vals}

        # (2) fetch witness limbs
        w_limbs_idx = witness_map[var_idx]
        print("w_limbs_idx", w_limbs_idx)

        # (3) multiply ⇒ term limbs
        term_limbs = add_multilimb_multiplication(builder, m_idx, c_limbs, w_limbs_idx)

        # (4) accumulate
        result_limbs.append(term_limbs)

    return result_limbs


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
        # add_range_check_constraints(builder, limb_idx)

        witness_map[j] = limb_idx

    # ------------------------------------------------------------------
    # 3. Constraint transformation – row by row
    # ------------------------------------------------------------------
    num_constraints = len(A)
    assert num_constraints == len(B) == len(C)

    for k in range(num_constraints):
        # Decompose *coefficients* ≥ bound (they are constants, so we keep them
        # as values for now).  Coefficients < bound stay single limb.
        rowA = A[k]
        rowB = B[k]
        rowC = C[k]

        a_k_limbs = {
            j: ([coeff] if coeff < norm_bound else decompose(coeff, norm_bound))
            for j, coeff in enumerate(rowA)
            if coeff
        }
        print("a_k_limbs", a_k_limbs)
        b_k_limbs = {
            j: ([coeff] if coeff < norm_bound else decompose(coeff, norm_bound))
            for j, coeff in enumerate(rowB)
            if coeff
        }
        print("b_k_limbs", b_k_limbs)
        c_k_limbs = {
            j: ([coeff] if coeff < norm_bound else decompose(coeff, norm_bound))
            for j, coeff in enumerate(rowC)
            if coeff
        }
        print("c_k_limbs", c_k_limbs)

        # Build limb representation of the dot products
        A_limbs = add_multilimb_dot_product(builder, 0, a_k_limbs, witness_map)
        B_limbs = add_multilimb_dot_product(builder, 1, b_k_limbs, witness_map)
        C_limbs = add_multilimb_dot_product(builder, 2, c_k_limbs, witness_map)

        # At this point we have the limbs of the dot products for R1CS row k of each matrix.
        # Now we need to add the constraints that enforce the Hadamard product:
        # (A·w) ⊙ (B·w) = C·w for the row k.
        add_multilimb_hadamard_constraint(builder, A_limbs, B_limbs, C_limbs)

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


def choose_decomposition_base(A, B, C, witness, q: int, norm_bound: int) -> int:
    """Return the largest *power-of-two* base ``d`` (< ``norm_bound``).

    Assumptions
    1. ``norm_bound`` itself is a power of two.
    2. We only decompose numbers in base ``d = 2^k``.

    The procedure below is still *conservative*: it checks that **all**
    powers ``d^e  (mod q)`` that can surface in the transformed system stay
    below ``norm_bound``.  The exponent bound is

        max_exp = 2 * max_bitlen - 1

    which safely covers dot-products and products of two numbers whose
    limbs were cut by ``norm_bound``.
    """

    assert norm_bound & (norm_bound - 1) == 0, "norm_bound must be a power of two"

    # --------------------------------------------------------------
    # 1. Maximum bit-length among all inputs (same as before)
    # --------------------------------------------------------------

    def _iter_values(mat):
        for row in mat:
            if isinstance(row, dict):
                yield from row.values()
            else:
                yield from row

    all_vals = (
        list(_iter_values(A))
        + list(_iter_values(B))
        + list(_iter_values(C))
        + list(witness)
    )

    max_val = max(abs(int(v)) for v in all_vals) if all_vals else 0

    max_bitlen = max_val.bit_length() if max_val else 1

    # --------------------------------------------------------------
    # 2. Iterate over powers of two < norm_bound, descending
    # --------------------------------------------------------------

    candidate = norm_bound >> 1  # highest power of two strictly < norm_bound

    while candidate >= 2:
        # Candidate-specific exponent bound ---------------------------
        log2d = int(math.log2(candidate))  # since candidate is power-of-two
        digits_needed = math.ceil(max_bitlen / log2d)
        max_exp_d = 2 * digits_needed - 1

        ok = True
        val = 1
        for _ in range(max_exp_d):
            val = (val * candidate) % q
            if val >= norm_bound:
                ok = False
                break

        if ok:
            return candidate

        candidate >>= 1  # next lower power of two

    raise ValueError(
        "Could not find a suitable power-of-two base – norm_bound too small relative to input values"
    )


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


def _sparse_rows_to_dense(rows: list[dict[int, int]], num_cols: int) -> list[list[int]]:
    dense_rows: list[list[int]] = []
    for row in rows:
        dense = [0] * num_cols
        for j, coeff in row.items():
            dense[j] = int(coeff)
        dense_rows.append(dense)

    return dense_rows


def print_dense_r1cs(
    A: list[dict[int, int]],
    B: list[dict[int, int]],
    C: list[dict[int, int]],
    witness: list[int],
) -> None:
    num_cols = len(witness)

    A_dense = _sparse_rows_to_dense(A, num_cols)
    B_dense = _sparse_rows_to_dense(B, num_cols)
    C_dense = _sparse_rows_to_dense(C, num_cols)

    print("Dense A matrix:")
    for row in A_dense:
        print(row)

    print("\nDense B matrix:")
    for row in B_dense:
        print(row)

    print("\nDense C matrix:")
    for row in C_dense:
        print(row)

    print("\nDense witness:")
    print(witness)


if __name__ == "__main__":

    q = 2**61 - 1
    norm_bound = 2**54

    # with open("mimc_r1cs.json", "r") as f:
    #     data = json.load(f)
    #     A = [[int(x, 16) for x in row] for row in data["A"]]
    #     B = [[int(x, 16) for x in row] for row in data["B"]]
    #     C = [[int(x, 16) for x in row] for row in data["C"]]
    #     witness = [int(x, 16) for x in data["witness"]]
    #
    # verify_r1cs_constraints(A, B, C, witness, q)
    # A_p, B_p, C_p, w_p = transform_r1cs_to_low_norm(A, B, C, witness, q, norm_bound)
    # verify_r1cs_constraints(A_p, B_p, C_p, w_p, q)
    #
    # # Serialize the transformed R1CS to JSON
    # transformed_data = {
    #     "A": [[hex(x) for x in row] for row in A_p],
    #     "B": [[hex(x) for x in row] for row in B_p],
    #     "C": [[hex(x) for x in row] for row in C_p],
    #     "witness": [hex(x) for x in w_p],
    # }
    #
    # with open("mimc_r1cs_transformed.json", "w") as f:
    #     json.dump(transformed_data, f, indent=2)
    # print("All constraints satisfied!")
    # print("#constraints:", len(A_p))
    # print("#witness variables:", len(w_p))
    # simple multiplication test

    witness = [1, 2**55, 3, 3 * 5 * 7 * 2**55]
    A = [[0, 5, 0, 0]]
    B = [[0, 0, 7, 0]]
    C = [[0, 0, 0, 1]]
    # base = choose_decomposition_base(A, B, C, witness, q, norm_bound)
    # print(f"base: {base}, 2**{int(math.log2(base))}")
    base = norm_bound

    verify_r1cs_constraints(A, B, C, witness, q)
    A_p, B_p, C_p, w_p = transform_r1cs_to_low_norm(A, B, C, witness, q, base)
    print("w_p:", w_p)
    print("A_p:", A_p)
    print("B_p:", B_p)
    print("C_p:", C_p)
    print_dense_r1cs(A_p, B_p, C_p, w_p)
    assert all(x < norm_bound for x in w_p)
    assert all(coeff <= norm_bound for row in A_p for coeff in row.values())
    assert all(coeff <= norm_bound for row in B_p for coeff in row.values())
    assert all(coeff <= norm_bound for row in C_p for coeff in row.values())
    verify_r1cs_constraints(A_p, B_p, C_p, w_p, q)

    print("All constraints satisfied!")
    print("#constraints:", len(A_p))
    print("#witness variables:", len(w_p))

    # ------------------------------------------------------------------
    # Additional test 1 – large witness values
    # ------------------------------------------------------------------

    a = 2**60 - 123  # ≥ norm_bound, triggers decomposition
    b = 2**60 - 45678  # ≥ norm_bound, triggers decomposition
    prod = (a * b) % q

    witness2 = [1, a, b, prod]
    A2 = [[0, 1, 0, 0]]
    B2 = [[0, 0, 1, 0]]
    C2 = [[0, 0, 0, 1]]

    verify_r1cs_constraints(A2, B2, C2, witness2, q)
    A2_p, B2_p, C2_p, w2_p = transform_r1cs_to_low_norm(
        A2, B2, C2, witness2, q, norm_bound
    )
    print_dense_r1cs(A2_p, B2_p, C2_p, w2_p)
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

    witness3 = [1, x, y]

    # A row applies big_coeff to x, B row multiplies by constant wire 1, C row expects the product
    A3 = [[0, big_coeff, 0]]  # coef on x
    B3 = [[1, 0, 0]]  # constant 1 (wire 0)
    C3 = [[0, 0, 1]]  # coef on y

    verify_r1cs_constraints(A3, B3, C3, witness3, q)
    A3_p, B3_p, C3_p, w3_p = transform_r1cs_to_low_norm(
        A3, B3, C3, witness3, q, norm_bound
    )
    verify_r1cs_constraints(A3_p, B3_p, C3_p, w3_p, q)

    print("\nLarge coefficient test passed!")
    print("#constraints:", len(A3_p))
    print("#witness variables:", len(w3_p))
