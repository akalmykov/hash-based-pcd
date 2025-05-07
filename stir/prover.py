import json
import math
import galois
from constants import FIELD192
from domain import Domain
import time
from fiat_schamir import Blake3Sponge
from folding import poly_fold
from pow import proof_of_work
from utils import stack_evals
from merkle import MerkleTree


class Prover:

    def calc_repetitions(self, log_inv_rate):
        # Using Provable soundness type with constant = 2
        return math.ceil(
            (self.soundness_type * self.protocol_security_level) / log_inv_rate
        )

    def pow_bits(self, log_inv_rate):
        repetitions = self.calc_repetitions(log_inv_rate)
        print(f"*** ??? log_inv_rate = {log_inv_rate}")
        print(f"*** ??? repetitions = {repetitions}")
        # TODO: This will change with eta
        scaling_factor = self.soundness_type
        achieved_security_bits = (log_inv_rate // scaling_factor) * repetitions
        remaining_security_bits = self.security_level - achieved_security_bits

        if remaining_security_bits <= 0:
            return 0
        else:
            return math.ceil(remaining_security_bits)  # Equivalent to ceil

    def __init__(self):
        start_time = time.time()
        self.field = galois.GF(FIELD192)
        end_time = time.time()
        print(f"Field init: {end_time - start_time} seconds")

        self.starting_degree = 2**10  # Evaluate polynomial: 296.5279061794281 seconds
        self.stopping_degree = 2**2
        self.field_size = 191
        self.protocol_security_level = 106
        self.security_level = 128
        self.starting_rate = 2
        self.folding_factor = 16
        self.ood_samples = 2
        self.soundness_type = 1  # Provable = 2, Conjecture = 1
        self.sponge = Blake3Sponge(self.field, 191)

        self.number_of_rounds, self.degrees = self.num_rounds_and_degrees()

        self.rates = [self.starting_rate]
        self.log_folding = self.folding_factor.bit_length() - 1
        self.rates.extend(
            [
                self.starting_rate + x * (self.log_folding - 1)
                for x in range(1, self.number_of_rounds + 1)
            ]
        )
        self.pow_bits = [self.pow_bits(x) for x in self.rates]
        self.repetitions = [self.calc_repetitions(x) for x in self.rates]

        print("self.folding_factor", self.folding_factor)
        # Note, this skips the last repetition
        for i in range(self.number_of_rounds):
            print("self.degrees[i]", self.degrees[i])
            self.repetitions[i] = min(
                self.repetitions[i], self.degrees[i] / self.folding_factor
            )

        assert self.number_of_rounds + 1 == len(self.rates)
        assert self.number_of_rounds + 1 == len(self.repetitions)

        print("Prover settings:")
        print(f"field: {self.field}")
        print(f"starting_degree: {self.starting_degree}")
        print(f"stopping_degree: {self.stopping_degree}")
        print(f"field_size: {self.field_size}")
        print(f"protocol_security_level: {self.protocol_security_level}")
        print(f"security_level: {self.security_level}")
        print(f"starting_rate: {self.starting_rate}")
        print(f"folding_factor: {self.folding_factor}")
        print(f"ood_samples: {self.ood_samples}")
        print(f"soundness_type: {self.soundness_type}")
        print(f"number_of_rounds: {self.number_of_rounds}")
        print(f"degrees: {self.degrees}")
        print(f"rates: {self.rates}")
        print(f"log_folding: {self.log_folding}")
        print(f"pow_bits: {self.pow_bits}")
        print(f"repetitions: {self.repetitions}")

    def num_rounds_and_degrees(self):
        d = self.starting_degree
        num_rounds = 0
        degrees = [d]
        while d > self.stopping_degree:
            assert (
                d % self.folding_factor == 0
            ), f"Degree {d} is not divisible by folding factor {self.folding_factor}"
            d //= self.folding_factor
            degrees.append(d)
            num_rounds += 1

        num_rounds -= 1
        degrees.pop()
        return num_rounds, degrees

    def commit(self, witness_polynomial):
        self.witness_polynomial = witness_polynomial
        start_time = time.time()
        self.domain = Domain(self.field, self.starting_degree, self.starting_rate)
        end_time = time.time()
        print(f"Domain: {end_time - start_time} seconds")

        start_time = time.time()
        evals = self.domain.evaluate_poly_coeff(self.witness_polynomial)
        end_time = time.time()
        print(f"Evaluate polynomial: {end_time - start_time} seconds")

        start_time = time.time()
        self.folded_evals = stack_evals(evals, self.folding_factor)
        end_time = time.time()
        print(f"Stack evals: {end_time - start_time} seconds")

        self.merkle_tree = MerkleTree(self.folded_evals)

    def round(self, folding_randomness, round_number):
        g_poly = poly_fold(
            self.witness_polynomial, self.folding_factor, folding_randomness
        )
        old_size = self.domain.get_size()
        self.domain.scale_with_offset(2)
        g_eval = self.domain.evaluate_poly_coeff(g_poly)
        g_folded_eval = stack_evals(g_eval, self.folding_factor)
        g_merkle = MerkleTree(g_folded_eval)
        self.sponge.absorb(bytes(g_merkle.root()))

        ood_rand = self.sponge.squeeze_several_f(self.ood_samples)

        _g_poly_evaluator = galois.Poly(list(reversed(g_poly)), field=self.field)
        betas = _g_poly_evaluator(ood_rand)
        self.sponge.absorb_field_vec(betas)
        comb_randomness = self.sponge.squeeze_f()
        next_folding_randomness = self.sponge.squeeze_f()
        print("next_folding_randomness", next_folding_randomness)
        scaling_factor = old_size // self.folding_factor
        num_repetitions = self.repetitions[round_number]
        stir_randomness_indexes = sorted(
            list(
                set(
                    self.sponge.squeeze_int(scaling_factor)
                    for i in range(num_repetitions)
                )
            )
        )
        pow_nonce = proof_of_work(self.sponge, self.pow_bits[round_number])
        _ = self.sponge.squeeze_f()
        queries_to_prev_ans = [self.folded_evals[i] for i in stir_randomness_indexes]
        print("scaling_factor", scaling_factor)
        queries_to_prev_proof = self.merkle_tree.generate_multi_proof(
            stir_randomness_indexes
        )
        queries_to_prev = (queries_to_prev_ans, queries_to_prev_proof)
        scaled_domain = self.domain.new_scaled(self.folding_factor)
        stir_randomness = [scaled_domain.element(i) for i in stir_randomness_indexes]
        # print("stir_randomness", stir_randomness)
        beta_answers = list(zip(ood_rand, betas))
        quotient_answers = [(x, _g_poly_evaluator(x)) for x in stir_randomness]
        quotient_answers.extend(beta_answers)

        quotient_set = list(ood_rand) + stir_randomness

        x_values, y_values = zip(*quotient_answers)
        # Convert to Galois Arrays
        x_values = self.field(list(x_values))
        y_values = self.field(list(y_values))

        ans_polynomial = galois.lagrange_poly(x_values, y_values)
        shake_polynomial = galois.Poly.Zero(field=self.field)
        for x, y in quotient_answers:
            num_polynomial = ans_polynomial - galois.Poly([y], field=self.field)
            den_polynomial = galois.Poly([1, -x], field=self.field)
            shake_polynomial = shake_polynomial + num_polynomial // den_polynomial

        vanishing_polynomial = galois.Poly.One(self.field)
        for s in quotient_set:
            vanishing_polynomial *= x - s
        numerator = _g_poly_evaluator + ans_polynomial
        quotient_polynomial = numerator // vanishing_polynomial
        print("quotient_polynomial", quotient_polynomial)

    def prove(self):
        self.sponge.absorb(bytes(self.merkle_tree.root()))
        folding_randomness = self.sponge.squeeze_f()

        proofs = []
        for round_number in range(self.number_of_rounds):
            round_proof = self.round(folding_randomness, round_number)
            proofs.append(round_proof)


if __name__ == "__main__":
    prover = Prover()
    with open("../../stir/poly_coeffs.json", "r") as f:
        fixed_poly = json.load(f)
    poly_coeffs = prover.field([int(fixed_poly[i]) for i in range(len(fixed_poly))])

    prover.commit(poly_coeffs)
    prover.prove()
