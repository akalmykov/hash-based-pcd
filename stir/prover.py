import json
import galois
from constants import FIELD192
from domain import Domain
import time
from fiat_schamir import Blake3Sponge
from folding import poly_fold
from utils import stack_evals
from merkle import MerkleTree


class Prover:
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
        self.sponge = Blake3Sponge(self.field, 191)

    def num_rounds(self):
        d = self.starting_degree
        num_rounds = 0
        while d > self.stopping_degree:
            assert (
                d % self.folding_factor == 0
            ), f"Degree {d} is not divisible by folding factor {self.folding_factor}"
            d //= self.folding_factor
            num_rounds += 1

        num_rounds -= 1
        return num_rounds

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

    def round(self, folding_randomness):
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

        betas = galois.Poly(list(reversed(g_poly)), field=self.field)(ood_rand)
        self.sponge.absorb_field_vec(betas)
        comb_randomness = self.sponge.squeeze_f()
        next_folding_randomness = self.sponge.squeeze_f()
        scaling_factor = old_size // self.folding_factor

    def prove(self):
        self.sponge.absorb(bytes(self.merkle_tree.root()))
        folding_randomness = self.sponge.squeeze_f()

        n = self.num_rounds()
        proofs = []
        for _ in range(n):
            round_proof = self.round(folding_randomness)
            proofs.append(round_proof)


if __name__ == "__main__":
    prover = Prover()
    with open("../../stir/poly_coeffs.json", "r") as f:
        fixed_poly = json.load(f)
    poly_coeffs = prover.field([int(fixed_poly[i]) for i in range(len(fixed_poly))])

    prover.commit(poly_coeffs)
    prover.prove()
