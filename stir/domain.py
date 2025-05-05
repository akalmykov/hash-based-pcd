import json
import galois
from constants import FIELD192, TEST_FIELD
import numpy as np
from utils import stack_evals, is_power_of_two, next_power_of_two, get_two_adicity
from merkle import MerkleTree, hash_field_vector, hash_pair, serialize_field_vector


# A function for extending a short array to a specified length:
def extendarray(short_array, length):
    size = len(short_array)
    # Validating short array is not too big:
    assert size <= length
    long_array = [0] * length
    # Write short_array contents into beginning of long_array:
    long_array[:size] = short_array[:size]
    return long_array


class Radix2EvaluationDomain:
    def __init__(self, GF, num_coeffs):
        self.GF = GF
        self.size = (
            num_coeffs if is_power_of_two(num_coeffs) else next_power_of_two(num_coeffs)
        )
        # print("self.size", self.size)
        self.log_size = self.size.bit_length() - 1  # size = 2 ^ log_size
        # print("self.log_size", self.log_size)
        if self.log_size > get_two_adicity(GF):
            raise ValueError("log_size is too large")
        # Compute the generator for the multiplicative subgroup.
        # It should be the 2^(log_size_of_group) root of unity.

        self.group_gen = GF.primitive_root_of_unity(self.size)
        self.size_as_field_element = GF(self.size)
        self.size_inv = self.size_as_field_element**-1
        self.group_gen_inv = self.group_gen**-1
        self.offset = GF(1)
        self.offset_inv = GF(1)
        self.offset_pow_size = GF(1)
        self.vandermonde_matrix = None
        # self.vandermonde_matrix = GF.Vandermonde(self.group_gen, self.size, self.size)
        # self.omega_pows = self.GF([self.group_gen ** (i) for i in range(self.size)])

    def omega_pows(self):
        return self.GF([self.group_gen ** (i) for i in range(self.size)])

    # # Method to construct evaluation matrix for the given domain
    # def as_vandermonde_matrix(self):
    #     self.evaluationmatrix = self.GF(
    #         [
    #             [self.group_gen ** (i * j) for i in range(self.size)]
    #             for j in range(self.size)
    #         ]
    #     )
    #     return self.evaluationmatrix

    # Method to evaluate a given polynomial over the given domain
    def evaluate_poly(self, poly):
        return self.evaluate_poly_coeff(list(reversed(poly.coeffs)))

    def evaluate_poly_coeff(self, poly_coeff):
        print(f"evaluating poly coeff size {len(poly_coeff)}, domain size {self.size}")
        return galois.ntt(poly_coeff, self.size, self.GF.order)


class Domain:

    def __init__(self, GF, degree, log_rho_inv):
        self.GF = GF
        self.degree = degree
        self.log_rho_inv = log_rho_inv
        self.size = degree * (1 << log_rho_inv)
        self.domain = Radix2EvaluationDomain(GF, self.size)

    def evaluation_domain_size(self):
        return self.domain.size

    def evaluate_poly(self, poly):
        return self.domain.evaluate_poly(poly)

    def evaluate_poly_coeff(self, poly_coeff):
        return self.domain.evaluate_poly_coeff(poly_coeff)

    def as_vandermonde_matrix(self):
        return self.domain.vandermonde_matrix

    def gen(self):
        return self.domain.group_gen

    def omega_pows(self):
        return self.domain.omega_pows()


if __name__ == "__main__":
    #   STIR
    # Targeting 128-bits of security - protocol running at 106-bits - soundness: Conjecture
    # Starting degree: 2^20, stopping_degree: 2^6
    # Starting rate: 2^-2, folding_factor: 16
    # Number of rounds: 3. OOD samples: 2
    # Rates: 2^-2, 2^-5, 2^-8, 2^-11
    # PoW bits: [22, 18, 16, 18]
    # Repetitions: [53, 22, 14, 10]
    # domain<f> size: 4194304
    # backing_domain: Radix2(Radix-2 multiplicative subgroup of size 4194304)
    # backing_domain.size(): 4194304
    # backing_domain.log_size_of_group(): 22

    # GF = galois.GF(FIELD192)
    # print(get_two_adicity(GF))
    # starting_rate = 2
    # starting_degree = 2**20
    # domain = Domain(GF, starting_degree, starting_rate)
    # print(domain.size)

    GF = galois.GF(17)
    print(GF.properties)
    domain = Domain(GF, 3, 2)
    domain_matrix = domain.as_vandermonde_matrix()
    print(domain_matrix)
    print(domain.gen())
    matrix = GF(
        [
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 3, 9, 10, 13, 5, 15, 11, 16, 14, 8, 7, 4, 12, 2, 6],
            [1, 9, 13, 15, 16, 8, 4, 2, 1, 9, 13, 15, 16, 8, 4, 2],
            [1, 10, 15, 14, 4, 6, 9, 5, 16, 7, 2, 3, 13, 11, 8, 12],
            [1, 13, 16, 4, 1, 13, 16, 4, 1, 13, 16, 4, 1, 13, 16, 4],
            [1, 5, 8, 6, 13, 14, 2, 10, 16, 12, 9, 11, 4, 3, 15, 7],
            [1, 15, 4, 9, 16, 2, 13, 8, 1, 15, 4, 9, 16, 2, 13, 8],
            [1, 11, 2, 5, 4, 10, 8, 3, 16, 6, 15, 12, 13, 7, 9, 14],
            [1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16, 1, 16],
            [1, 14, 9, 7, 13, 12, 15, 6, 16, 3, 8, 10, 4, 5, 2, 11],
            [1, 8, 13, 2, 16, 9, 4, 15, 1, 8, 13, 2, 16, 9, 4, 15],
            [1, 7, 15, 3, 4, 11, 9, 12, 16, 10, 2, 14, 13, 6, 8, 5],
            [1, 4, 16, 13, 1, 4, 16, 13, 1, 4, 16, 13, 1, 4, 16, 13],
            [1, 12, 8, 11, 13, 3, 2, 7, 16, 5, 9, 6, 4, 14, 15, 10],
            [1, 2, 4, 8, 16, 15, 13, 9, 1, 2, 4, 8, 16, 15, 13, 9],
            [1, 6, 2, 12, 4, 7, 8, 14, 16, 11, 15, 5, 13, 10, 9, 3],
        ]
    )
    # assert (domain_matrix == matrix).all()
    p = galois.Poly([2, 7, 5, 3], field=GF)
    print(p)
    print(p(domain.omega_pows()))
    print(domain.evaluate_poly(p))
    assert (p(domain.omega_pows()) == domain.evaluate_poly(p)).all()
    print(domain.size)
    print(domain.evaluation_domain_size())
    print(galois.ntt(p.coeffs, GF.order - 1, GF.order))
    print(galois.ntt(list(reversed([2, 7, 5, 3])), GF.order - 1, GF.order))
    print(galois.ntt(list(reversed(p.coeffs)), GF.order - 1, GF.order))
    print(np.fft.fft(p.coeffs, n=GF.order - 1))

    with open("../../stir/domain_points.json", "r") as f:
        domain_points = json.load(f)

    GF = galois.GF(FIELD192)
    t = GF([int(domain_points[i]) for i in range(len(domain_points))])

    print(GF.properties)
    domain = Domain(GF, 2**10, 2)
    print(len(domain_points))
    print(len(domain.omega_pows()))
    assert len(domain_points) == len(domain.omega_pows())
    assert (t == domain.omega_pows()).all()

    with open("../../stir/evals.json", "r") as f:
        evals = json.load(f)
    t_evals = GF([int(evals[i]) for i in range(len(evals))])
    print(t_evals[-1])

    with open("../../stir/poly_coeffs.json", "r") as f:
        poly_coeffs = json.load(f)
    t_poly_coeffs = GF([int(poly_coeffs[i]) for i in range(len(poly_coeffs))])

    my_evals = domain.evaluate_poly_coeff(t_poly_coeffs)
    assert (my_evals == t_evals).all()
    # assert (domain.omega_pows() == domain_points).all()
    stacked_evals = stack_evals(my_evals, 16)

    with open("../../stir/folded_evals.json", "r") as f:
        folded_evals_json = json.load(f)

    t_folded_evals = []
    for leaf in folded_evals_json:
        t_folded_evals.append(GF([GF(leaf[i]) for i in range(len(leaf))]))

    for i in range(len(stacked_evals)):
        assert (stacked_evals[i] == t_folded_evals[i]).all()

    with open("../../stir/merkle_tree.json", "r") as f:
        merkle_tree = json.load(f)

    leaf_nodes = []
    for node in merkle_tree["leaf_nodes"]:
        # Extract the byte array from the string format "SHA3Digest([...])"
        byte_str = node[11:-1]  # Remove "SHA3Digest(" and ")"
        bytes_list = [int(x) for x in byte_str.strip("[]").split(", ")]
        leaf_nodes.append(bytes_list)
    from hashlib import sha3_256

    last_eval = t_folded_evals[-1]
    last_eval_bytes = serialize_field_vector(last_eval)

    # Calculate SHA3-256 digest
    sha3 = sha3_256()
    sha3.update(last_eval_bytes)
    calculated_digest = list(sha3.digest())
    # print(leaf_nodes)
    # print(calculated_digest)
    assert leaf_nodes[-1] == calculated_digest
    assert calculated_digest == hash_field_vector(last_eval)
    assert len(leaf_nodes) == len(stacked_evals)
    for i, leaf in enumerate(leaf_nodes):
        assert leaf == hash_field_vector(stacked_evals[i])
    # print("Bytes in last_eval_bytes:", [b for b in last_eval_bytes])

    # l1 = [2568590649130636328156454025689191507038995861110281789944]
    # t1 = b"".join([int(x).to_bytes(24, "little") for x in l1])
    # print("Bytes in t1:", [b for b in t1])

    # calculate non leaf nodes from staked_evals
    # store the non-leaf nodes in level order. The first element is the root node.
    # the ith nodes (starting at 1st) children are at indices 2*i, 2*i+1
    non_leaf_nodes = []
    for node in merkle_tree["non_leaf_nodes"]:
        # Extract the byte array from the string format "SHA3Digest([...])"
        byte_str = node[11:-1]  # Remove "SHA3Digest(" and ")"
        bytes_list = [int(x) for x in byte_str.strip("[]").split(", ")]
        non_leaf_nodes.append(bytes_list)

    right = hash_field_vector(stacked_evals[-1])
    left = hash_field_vector(stacked_evals[-2])
    last_non_leaf = hash_pair(left, right)

    mt = MerkleTree(stacked_evals)
    for i, leaf in enumerate(mt.leaf_nodes):
        assert leaf == hash_field_vector(stacked_evals[i])

    print("mt.non_leaf_nodes[-1]", mt.non_leaf_nodes[-1])
    print("last_non_leaf", last_non_leaf)
    assert mt.non_leaf_nodes[-1] == last_non_leaf

    j = 0
    for i in range(0, len(stacked_evals), 2):
        h1 = hash_pair(
            hash_field_vector(stacked_evals[i]),
            hash_field_vector(stacked_evals[i + 1]),
        )
        h2 = mt.non_leaf_nodes[len(stacked_evals) // 2 - 1 + j]
        assert h1 == h2
        j += 1

    for i in range(len(mt.non_leaf_nodes)):
        assert mt.non_leaf_nodes[i] == non_leaf_nodes[i]
    print("Merkle tree root: ", non_leaf_nodes[0])
    # l2 = [2568590649130636328156454025689191507038995861110281789944]
    # t2 = b"".join([int(x).to_bytes(24, "big") for x in l2])
    # print("Bytes in t2:", [b for b in t2])

    # l3 = [2568590649130636328156454025689191507038995861110281789944]
    # t3 = b"".join([int(x).to_bytes(24, "big") for x in l3])
    # print("Bytes in t3:", [b for b in t3])

    pass
