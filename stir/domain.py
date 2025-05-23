import json
import galois
from .constants import FIELD192, TEST_FIELD
import numpy as np
from .utils import stack_evals, is_power_of_two, next_power_of_two, get_two_adicity
from .merkle import MerkleTree, hash_field_vector, hash_pair, serialize_field_vector
from dataclasses import dataclass


@dataclass
class EvaluationDomainConfig:
    size: int
    size_as_field_element: object
    size_inv: object
    group_gen: object
    group_gen_inv: object
    offset: object
    offset_inv: object
    offset_pow_size: object


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
    @classmethod
    def from_config(cls, GF, config):
        """
        Create a Radix2EvaluationDomain instance from a configuration object.

        Args:
            GF: The Galois Field
            config: EvaluationDomainConfig with domain parameters
        """
        instance = cls.__new__(cls)
        instance.GF = GF

        # Copy fields from config
        instance.size = config.size
        instance.size_as_field_element = config.size_as_field_element
        instance.size_inv = config.size_inv
        instance.group_gen = config.group_gen
        instance.group_gen_inv = config.group_gen_inv
        instance.offset = config.offset
        instance.offset_inv = config.offset_inv
        instance.offset_pow_size = config.offset_pow_size

        # Compute other fields as in original constructor
        instance.log_size = instance.size.bit_length() - 1
        if instance.log_size > get_two_adicity(GF):
            raise ValueError("log_size is too large")

        return instance

    def __init__(self, GF, num_coeffs) -> None:
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

    def omega_pows(self):
        return self.GF([self.offset * self.group_gen ** (i) for i in range(self.size)])

    # Method to evaluate a given polynomial over the given domain
    def evaluate_poly(self, poly):
        return self.evaluate_poly_coeff(list(reversed(poly.coeffs)))

    def coset_transform(self, coeffs):
        size = self.size
        coeffs = coeffs.copy()
        first_chunk = coeffs[:size].copy()

        offset = self.offset

        for i in range(1, len(coeffs) // size):
            chunk = coeffs[i * size : (i + 1) * size]
            if offset == 1:
                first_chunk += chunk
            else:
                offset_power = offset ** (i * size)
                first_chunk += offset_power * chunk
        return first_chunk

    def naive_eval(self, poly_coeff):
        pows = self.omega_pows()
        p = galois.Poly(list(reversed(poly_coeff)), field=self.GF)
        return p(pows)

    def evaluate_poly_coeff(self, poly_coeff):
        print(f"evaluating poly coeff size {len(poly_coeff)}, domain size {self.size}")
        if self.offset > 1:
            coeffs = []
            for i, coeff in enumerate(poly_coeff):
                coeffs.append(coeff * self.offset ** (i))
            return galois.ntt(
                coeffs,
                self.size,
                self.GF.order,
            )
        else:
            return galois.ntt(poly_coeff, self.size, self.GF.order)

    def scale(self, power):
        # TODO: refactor DRY
        self.size = self.size // power
        self.log_size = self.size.bit_length() - 1
        if self.log_size > get_two_adicity(self.GF):
            raise ValueError("log_size is too large")

        self.size_as_field_element = self.GF(self.size)
        self.size_inv = self.size_as_field_element**-1
        group_gen = self.group_gen**power
        group_gen_inv = self.group_gen_inv**power
        self.offset = (self.offset**power) * self.group_gen
        self.offset_inv = (self.offset_inv**power) * self.group_gen_inv
        self.group_gen = group_gen
        self.group_gen_inv = group_gen_inv

        print("self.offset", self.offset)
        print("group_gen", group_gen)
        print("group_gen_inv", group_gen_inv)
        print("offset", self.offset)
        print("offset_inv", self.offset_inv)

    def scale_generator_by(self, power):
        self.size = self.size // power
        self.log_size = self.size.bit_length() - 1
        self.size_as_field_element = self.GF(self.size)
        self.size_inv = self.size_as_field_element**-1
        group_gen = self.group_gen**power
        group_gen_inv = group_gen**-1
        self.offset = self.offset**power
        self.offset_inv = self.offset_inv**power
        self.group_gen = group_gen
        self.group_gen_inv = group_gen_inv
        self.offset_pow_size = self.offset**self.size

    def element(self, i):
        result = self.group_gen**i
        if self.offset != self.GF(1):
            result *= self.offset
        return result

    def __str__(self):
        return (
            f"Radix2EvaluationDomain(\n"
            f"  GF: {self.GF}\n"
            f"  size: {self.size}\n"
            f"  log_size: {self.log_size}\n"
            f"  group_gen: {self.group_gen}\n"
            f"  group_gen_inv: {self.group_gen_inv}\n"
            f"  size_as_field_element: {self.size_as_field_element}\n"
            f"  size_inv: {self.size_inv}\n"
            f"  offset: {self.offset}\n"
            f"  offset_inv: {self.offset_inv}\n"
            f"  offset_pow_size: {self.offset_pow_size}\n"
            f")"
        )


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

    # def as_vandermonde_matrix(self):
    #     return self.domain.vandermonde_matrix

    def gen(self):
        return self.domain.group_gen

    def omega_pows(self):
        return self.domain.omega_pows()

    def _scale_with_offset(self, power):
        self.domain.scale(power)

    def new_scale_with_offset(self, power):
        scaled_domain = Radix2EvaluationDomain(self.GF, self.size)
        scaled_domain.scale(power)
        return scaled_domain

    def new_scaled(self, power):
        scaled_domain = Radix2EvaluationDomain(self.GF, self.size)
        scaled_domain.scale_generator_by(power)
        return scaled_domain

    def get_size(self):
        return self.domain.size

    def element(self, i):
        return self.domain.element(i)


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
    # domain_matrix = domain.as_vandermonde_matrix()
    # print(domain_matrix)
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

    with open("./test/original-stir-traces/domain_points.json", "r") as f:
        domain_points = json.load(f)

    GF = galois.GF(FIELD192)
    t = GF([int(domain_points[i]) for i in range(len(domain_points))])

    print(GF.properties)
    domain = Domain(GF, 2**10, 2)
    print(len(domain_points))
    print(len(domain.omega_pows()))
    assert len(domain_points) == len(domain.omega_pows())
    assert (t == domain.omega_pows()).all()

    with open("./test/original-stir-traces/evals.json", "r") as f:
        evals = json.load(f)
    t_evals = GF([int(evals[i]) for i in range(len(evals))])
    print(t_evals[-1])

    with open("./test/original-stir-traces/poly_coeffs.json", "r") as f:
        poly_coeffs = json.load(f)
    t_poly_coeffs = GF([int(poly_coeffs[i]) for i in range(len(poly_coeffs))])

    my_evals = domain.evaluate_poly_coeff(t_poly_coeffs)
    assert (my_evals == t_evals).all()
    # assert (domain.omega_pows() == domain_points).all()
    stacked_evals = stack_evals(my_evals, 16)

    with open("./test/original-stir-traces/folded_evals.json", "r") as f:
        folded_evals_json = json.load(f)

    t_folded_evals = []
    for leaf in folded_evals_json:
        t_folded_evals.append(GF([GF(leaf[i]) for i in range(len(leaf))]))

    for i in range(len(stacked_evals)):
        assert (stacked_evals[i] == t_folded_evals[i]).all()

    with open("./test/original-stir-traces/merkle_tree.json", "r") as f:
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
