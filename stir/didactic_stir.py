from sympy import ntt, intt
import numpy as np
import galois

# Security Warning: The finite field arithmetic implemented in galois are not constant-time, but were instead designed for performance. As such, the library could be vulnerable to a side-channel timing attack. This library is not intended for production security, but instead for research & development, reverse engineering, cryptanalysis, experimentation, and general education.

# Preliminaries

# We use galois to define _finite_field_array. This functions like a numpy array, except
# that the array consists of elements of the finite field of size _field_size.
_field_size = 97
_finite_field_array = galois.GF(_field_size)
_primitive_root = _finite_field_array.primitive_element
_trace_length = 8
_reed_solomon_expansion_factor = 4
_block_length = _trace_length * _reed_solomon_expansion_factor
_generator = _primitive_root ** ((_field_size - 1) // _block_length)


# An EvaluationDomain is a collection of points on which to evaluate a polynomial
class EvaluationDomain:
    # Construct by directly providing the domain as an array
    def __init__(self, array):
        self.size = len(array)
        self.domain = array
        if self.size > 1:
            self.evaluation_domain_of_half_size = EvaluationDomain(
                [self.domain[2 * i] for i in range(self.size // 2)]
            )

    # Alternate constructor
    @classmethod
    def from_shift_generator_size(self, shift, generator, size):
        return EvaluationDomain(
            np.array(
                [(shift * mod_pow(generator, i)) % _field_size for i in range(size)]
            )
        )

    # Method to construct evaluation matrix for the given domain
    def as_vandermonde_matrix(self):
        self.evaluationmatrix = np.array(
            [
                [mod_pow(self.domain[j], i) for i in range(self.size)]
                for j in range(self.size)
            ]
        )
        return self.evaluationmatrix

    # Method to evaluate a given polynomial over the given domain
    def evaluate(self, coeffs):
        self.commitment = (
            self.as_vandermonde_matrix() @ extendarray(coeffs, self.size) % _field_size
        )
        return self.commitment


# A function that performs modular exponentiation
def mod_pow(base, exp):
    return pow(int(base), exp, _field_size)


# A function that returns the reciprocal modulo _field_size
def field_reciprocal(x):
    # By Fermat's Little Theorem, x^95 is the reciprocal of x (for non-zero x).
    x = _finite_field_array(x % _field_size)
    xinv = _finite_field_array(1) / x
    return xinv


# A function that returns an element-wise reciprocal of an array, modulo _field_size
def reciprocal_elementwise(array):
    return [field_reciprocal(x) for x in array]


# A function for extending a short array to a specified length:
def extendarray(short_array, length):
    size = len(short_array)
    # Validating short array is not too big:
    assert size <= length
    long_array = [0] * length
    # Write short_array contents into beginning of long_array:
    long_array[:size] = short_array[:size]
    return long_array


# A function that evaluates a polynomial at a given input
def evaluate_at_point(coeffs, x):
    # We evaluate f(x) as the dot product of [coefficients] and [powers]=[1, x, x^2, ...]
    powers = np.array([mod_pow(x, i) for i in range(len(coeffs))])
    eval = np.dot(coeffs, powers)
    return eval


def main():
    print(f"Generator {_generator}")
    print(f"Primitive root {_primitive_root}")

    # Trace: Data Columns
    d1_trace = np.array([24, 30, 54, 84, 78, 15, 29, 50])
    print("Trace Data Columns \n", d1_trace)

    d1_coeffs = np.array(intt(d1_trace, prime=_field_size))

    print("Coefficients of Trace Polynomials: \n", d1_coeffs)
    folding_factor = 2
    bs08_matrix = to_coefficient_matrix(
        d1_coeffs, len(d1_coeffs) // folding_factor, folding_factor
    )
    print("BS08 Matrix: \n", bs08_matrix.matrix)
    folding_randomness = 5
    folded_polynomial = bs08_matrix.fold_by_col(folding_randomness)
    print("Folded Polynomial: \n", folded_polynomial)


if __name__ == "__main__":
    main()
