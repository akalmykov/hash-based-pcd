import galois
from constants import FIELD192, TEST_FIELD
import numpy as np


def is_power_of_two(n):
    return n & (n - 1) == 0


def next_power_of_two(n):
    return 1 << (n - 1).bit_length()


def get_two_adicity(GF):
    n = GF.order - 1

    # Count trailing zeros to get two-adicity
    adicity = 0
    while n & 1 == 0:
        adicity += 1
        n >>= 1

    return adicity


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
        return self.domain.omega_pows


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
    print(domain.evaluate(p))
    assert (p(domain.omega_pows()) == domain.evaluate(p)).all()
    print(domain.size)
    print(domain.evaluation_domain_size())
    print(galois.ntt(p.coeffs, GF.order - 1, GF.order))
    print(galois.ntt(list(reversed([2, 7, 5, 3])), GF.order - 1, GF.order))
    print(galois.ntt(list(reversed(p.coeffs)), GF.order - 1, GF.order))
    print(np.fft.fft(p.coeffs, n=GF.order - 1))

    pass
