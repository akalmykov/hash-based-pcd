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


# Define mix function for inside FRI
def mix(array, parameter):
    mixed_array = np.array([0] * len(array[0]))
    for i in range(len(array)):
        assert len(array[i]) == len(array[0])
        mixed_array = (
            mixed_array + (np.array(array[i]) * mod_pow(parameter, i))
        ) % _field_size
    return mixed_array


# Define quotient function for DEEP
def DEEPquotient(block, domain, testpoint, eval_testpoint):
    eval_quotient = (
        (block - eval_testpoint)
        * reciprocal_elementwise(domain - testpoint)
        % _field_size
    )
    return eval_quotient


# Define a function that splits a block into 2 blocks of half the length
def split(block, large_domain):
    small_domain = large_domain.evaluation_domain_of_half_size
    # Re-write block as a galois array
    block = _finite_field_array(block)
    evalmatrix_inv = np.linalg.inv(
        _finite_field_array(large_domain.as_vandermonde_matrix())
    )

    # Given a block u, we construct smaller blocks v_0, v_1

    # The defining relationship between u and v_i is easy to state in coefficient representation:
    #   Writing u_bar as the coefficient form of u, and v_bar_i as the coefficient form of v_i, we write:
    #   u_bar(x) = v_bar_0(x^2) + (x * v_bar_1(x^2))
    # We compute the coefficients of v_bar_0 and v_bar_1 through the following matrix multiplication:
    # _finite_field_array is the finite field of order _field_size:

    # Computing coefficient form of input:
    u_bar = evalmatrix_inv @ block

    # Writing coefficients of output polynomials:
    v_bar_0 = extendarray(
        [u_bar[2 * i] for i in range(len(u_bar) // 2)], len(small_domain.domain)
    )
    v_bar_1 = extendarray(
        [u_bar[1 + 2 * i] for i in range(len(u_bar) // 2)], len(small_domain.domain)
    )

    # Computing output blocks
    v_0 = small_domain.evaluate(v_bar_0)
    v_1 = small_domain.evaluate(v_bar_1)
    return [v_0, v_1]


def fri(block, input_domain, mixing_parameters, total_rounds, rounds_remaining):
    while rounds_remaining > 0:
        print("\n--- begin round", total_rounds - rounds_remaining + 1, "---")
        print("D_super_(", total_rounds - rounds_remaining, "):", input_domain.domain)
        split0, split1 = split(block, input_domain)
        next_domain = input_domain.evaluation_domain_of_half_size
        print("Split blocks: \n", split0, "\n", split1)
        print("Now mix those blocks using the mix parameter", mixing_parameters[0])
        round_commitment = mix([split0, split1], mixing_parameters[0])
        print(
            "FRI Commitment (aka f_super_(",
            total_rounds - rounds_remaining + 1,
            "):",
            round_commitment,
        )
        return fri(
            round_commitment,
            next_domain,
            mixing_parameters[1:],
            total_rounds,
            rounds_remaining - 1,
        )
    if rounds_remaining == 0:
        print("D_super_(", total_rounds - rounds_remaining, "):", input_domain.domain)
        return block


# Main ZKP Function:
def main():
    print("START \n")

    # A Didactic STARK
    # Proving computational integrity of a Fibonacci Program mod _field_size

    # We input 24 and 30 into a simple fibonacci program.
    # The inputs are loaded into d1_trace[0] and d2_trace[0]
    # We compute 4 rounds of fibonacci addition.
    # At each round, we write d3_trace[i] = d1_trace[i] + d2_trace[i]
    # After 4 rounds d3_trace[3] contains the output of the computation.

    # We claim that the following trace is a valid instantiation of this algorithm.
    # The trace is organized in columns, which are divided into two types: Data & Control.

    print(
        "This code is based on a Python implementation of the STARK by Hand Explainer."
    )
    print("https://dev.risczero.com/proof-system/stark-by-hand")
    print("")
    print("The Padded Execution Trace (See Lessons 1-3 in STARK by Hand)")
    print(f"Generator {_generator}")
    print(f"Primitive root {_primitive_root}")

    # Trace: Data Columns
    d1_trace = np.array([24, 30, 54, 84, 78, 15, 29, 50])
    d2_trace = np.array([30, 54, 84, 41, 2, 77, 21, 36])
    d3_trace = np.array([54, 84, 41, 28, 71, 17, 92, 33])
    print("Trace Data Columns \n", d1_trace, "\n", d2_trace, "\n", d3_trace, "\n")

    # Trace: Control Columns
    # Init steps are flagged in c1_trace
    # Computation steps are flagged in c2_trace
    # Termination step is flagged in c3_trace
    # 0s at the end of each control column correspond to the padding of the trace
    # c1_trace = np.array([1, 0, 0, 0, 0, 0, 0, 0])
    # c2_trace = np.array([0, 1, 1, 1, 0, 0, 0, 0])
    # c3_trace = np.array([0, 0, 0, 1, 0, 0, 0, 0])
    # print("Trace Control Columns \n", c1_trace, "\n", c2_trace, "\n", c3_trace, "\n")
    # trace_data = [d1_trace, d2_trace, d3_trace, c1_trace, c2_trace, c3_trace]

    # We will construct a zero-knowledge proof that:
    # this trace represents a program that satisfies these 6 rules:
    # 1) Fibonacci words here
    # 2) d1_trace[0] == 24  (init 1 constraint)
    # 3) d2_trace[0] == 30  (init 2 constraint)
    # 4) d3_trace[3] == 28  (termination constraint)
    # 5) if c2_trace[i] == 1, then d2_trace[i] == d1_trace[i+1]
    # 6) if c2_trace[i] == 1, then d3_trace[i] == d2_trace[i+1}

    print("Lesson 4: Constructing Trace Polynomials")
    # Trace Polynomials
    # We run an iNTT on each trace column to construct the associated trace polynomial
    d1_coeffs = np.array(intt(d1_trace, prime=_field_size))
    d2_coeffs = np.array(intt(d2_trace, prime=_field_size))
    d3_coeffs = np.array(intt(d3_trace, prime=_field_size))

    # c1_coeffs = np.array(intt(c1_trace, prime=_field_size))
    # c2_coeffs = np.array(intt(c2_trace, prime=_field_size))
    # c3_coeffs = np.array(intt(c3_trace, prime=_field_size))
    print(
        "Coefficients of Trace Polynomials 1,2,...,6: \n",
        d1_coeffs,
        "\n",
        d2_coeffs,
        "\n",
        d3_coeffs,
        "\n",
        # c1_coeffs,
        # "\n",
        # c2_coeffs,
        # "\n",
        # c3_coeffs,
        # "\n",
    )

    # Evaluating Trace Polynomials over the powers of 5^12 would return the original trace data
    # 5^12 is the first 8-th root of unity for F_97, sympy/intt uses the first root of unity
    # There are φ(8)=4 primitive 8-th roots of unity, where φ is Euler's totient function

    print("Verifying that polynomials are correct")
    expanded_domain = EvaluationDomain.from_shift_generator_size(1, 5**12, 8)
    print(expanded_domain.evaluate(d1_coeffs))
    print(expanded_domain.evaluate(d2_coeffs))
    print(expanded_domain.evaluate(d3_coeffs))
    print("Evaluate a linear combination of d1_coeffs + d2_coeffs")
    print(expanded_domain.evaluate(d1_coeffs + 2 * d2_coeffs))

    # Evaluating Trace Polynomials over the "expanded domain" gives a "trace block."
    expanded_domain = EvaluationDomain.from_shift_generator_size(1, 5**3, 32)
    d1_reedsolomonexpansion = expanded_domain.evaluate(d1_coeffs)
    d2_reedsolomonexpansion = expanded_domain.evaluate(d2_coeffs)
    d3_reedsolomonexpansion = expanded_domain.evaluate(d3_coeffs)
    # c1_reedsolomonexpansion = expanded_domain.evaluate(c1_coeffs)
    # c2_reedsolomonexpansion = expanded_domain.evaluate(c2_coeffs)
    # c3_reedsolomonexpansion = expanded_domain.evaluate(c3_coeffs)
    print(
        "Reed Solomon Expansion of Trace Blocks: \n",
        d1_reedsolomonexpansion,
        "\n",
        d2_reedsolomonexpansion,
        "\n",
        d3_reedsolomonexpansion,
        "\n",
        # c1_reedsolomonexpansion,
        # "\n",
        # c2_reedsolomonexpansion,
        # "\n",
        # c3_reedsolomonexpansion,
    )
    print("Evaluate an expanded linear combination of d1_coeffs and d2_coeffs")
    print(expanded_domain.evaluate(d1_coeffs + 2 * d2_coeffs))
    print(d1_reedsolomonexpansion + 2 * d2_reedsolomonexpansion)

    print("Note that every 4th entry matches the original trace data.")
    print("This is a degree 4 Reed Solomon expansion of the original trace. \n")
    trace_reedsolomonexpansion = np.array(
        [
            d1_reedsolomonexpansion,
            d2_reedsolomonexpansion,
            d3_reedsolomonexpansion,
            # c1_reedsolomonexpansion,
            # c2_reedsolomonexpansion,
            # c3_reedsolomonexpansion,
        ]
    )

    print("Lesson 5: ZK Commitments of the Trace Data")
    print(
        'To maintain a zero-knowledge protocol, the trace polynomials are evaluated over a "zk commitment domain", {5^1, 5^4, ..., 5^94}.'
    )
    # Evaluating Trace Polynomials over the zk commitment domain gives "zero-knowledge trace blocks."
    zk_commitment_domain = EvaluationDomain.from_shift_generator_size(5, 5**3, 32)
    d1_zkcommitment = zk_commitment_domain.evaluate(d1_coeffs)
    d2_zkcommitment = zk_commitment_domain.evaluate(d2_coeffs)
    d3_zkcommitment = zk_commitment_domain.evaluate(d3_coeffs)
    # c1_zkcommitment = zk_commitment_domain.evaluate(c1_coeffs)
    # c2_zkcommitment = zk_commitment_domain.evaluate(c2_coeffs)
    # c3_zkcommitment = zk_commitment_domain.evaluate(c3_coeffs)
    print(
        "Zero-Knowledge Trace Blocks: \n",
        d1_zkcommitment,
        "\n",
        d2_zkcommitment,
        "\n",
        d3_zkcommitment,
        "\n",
        # c1_zkcommitment,
        # "\n",
        # c2_zkcommitment,
        # "\n",
        # c3_zkcommitment,
    )
    print(
        "These zk-commitment blocks do not share any evaluation points with the original trace data."
    )
    print(
        "These evaluations are committed to Merkle Trees (one for the data blocks and one for the control blocks."
    )
    print(
        'The term "polynomial commitment" describes this process of a Merkle commitment where the leaves of the tree are evaluations of a polynomial. \n'
    )
    print(
        "The root of these two Merkle Trees are the first entries on the RISC Zero seal."
    )
    print(
        "The third entry on the RISC Zero seal is a `Accum` commitment; that commitment is not necessary in this simplified example. \n"
    )


if __name__ == "__main__":
    main()
