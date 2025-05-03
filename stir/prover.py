import galois
from constants import FIELD192
from domain import Domain
import time
from merkle import MerkleTree


def stack_evals(evals, folding_factor):
    assert len(evals) % folding_factor == 0
    size_of_the_new_domain = len(evals) // folding_factor
    stacked_evals = []
    for i in range(size_of_the_new_domain):
        new_evals = []
        for j in range(folding_factor):
            new_evals.append(evals[i + size_of_the_new_domain * j])
        stacked_evals.append(new_evals)
    return stacked_evals


class Prover:
    def __init__(self):
        start_time = time.time()
        self.field = galois.GF(FIELD192)
        end_time = time.time()
        print(f"Field init: {end_time - start_time} seconds")

        self.starting_degree = 2**20  # Evaluate polynomial: 296.5279061794281 seconds
        self.stopping_degree = 2**6
        self.protocol_security_level = 106
        self.security_level = 128
        self.starting_rate = 2
        self.folding_factor = 16

    def commit(self, witness_polynomial):
        start_time = time.time()
        domain = Domain(self.field, self.starting_degree, self.starting_rate)
        end_time = time.time()
        print(f"Domain: {end_time - start_time} seconds")

        start_time = time.time()
        evals = domain.evaluate_poly_coeff(witness_polynomial)
        end_time = time.time()
        print(f"Evaluate polynomial: {end_time - start_time} seconds")

        start_time = time.time()
        folded_evals = stack_evals(evals, self.folding_factor)
        end_time = time.time()
        print(f"Stack evals: {end_time - start_time} seconds")

        self.merkle_tree = MerkleTree(self.field, folded_evals)


if __name__ == "__main__":
    prover = Prover()
    poly = prover.field.Random(prover.starting_degree + 1, seed=1)
    prover.commit(poly)
