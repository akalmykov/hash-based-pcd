from .prover import Proof
from .stir_parameters import FullParameters
from .merkle import MerkleTree
from .fiat_shamir import Blake3Sponge
from .domain import Domain
from dataclasses import dataclass
from typing import Any
import galois
from enum import Enum
from typing import Tuple
from .prover import RoundProof
from .pow import proof_of_work_verify


class OracleType(Enum):
    INITIAL = "Initial"
    VIRTUAL = "Virtual"


@dataclass
class VerificationState:
    oracle: OracleType
    domain_gen: galois.FieldArray
    domain_size: int
    domain_offset: galois.FieldArray
    root_of_unity: galois.FieldArray
    folding_randomness: galois.FieldArray
    num_round: int


class Verifier:
    def __init__(self, GF, parameters: FullParameters):
        self.GF = GF
        self.parameters = parameters

    def verify(self, commitment: MerkleTree, proof: Proof):
        if len(proof.final_polynomial) > self.parameters.stopping_degree:
            return False
        current_root = commitment.root()
        for round_proof in proof.round_proofs:
            multipath = round_proof.queries_to_prev[1]
            if not multipath.verify(current_root, round_proof.queries_to_prev[0]):
                return False
            else:
                current_root = round_proof.g_root

        if not proof.queries_to_final[1].verify(
            current_root, proof.queries_to_final[0]
        ):
            return False

        sponge = Blake3Sponge(self.GF, 191)
        print("commitment.root()", commitment.root())
        sponge.absorb(bytes(commitment.root()))
        folding_randomness = sponge.squeeze_f()
        domain = Domain(
            self.GF, self.parameters.starting_degree, self.parameters.starting_rate
        )
        domain_gen = domain.element(1)
        domain_size = domain.get_size()
        print("folding_randomness:", folding_randomness)
        print("domain_gen:", domain_gen)
        print("domain_size:", domain_size)

        verification_state = VerificationState(
            oracle=OracleType.INITIAL,
            domain_gen=domain_gen,
            domain_size=domain_size,
            domain_offset=self.GF(1),
            root_of_unity=domain_gen,
            folding_randomness=folding_randomness,
            num_round=0,
        )

        for round_proof in proof.round_proofs:
            valid, new_state = self.round(sponge, round_proof, verification_state)
            if not valid:
                return False
            verification_state = new_state
        return True

    def round(
        self,
        sponge: Blake3Sponge,
        round_proof: RoundProof,
        verification_state: VerificationState,
    ) -> Tuple[bool, VerificationState]:
        sponge.absorb(round_proof.g_root)
        ood_randomness = sponge.squeeze_several_f(self.parameters.ood_samples)
        sponge.absorb_field_vec(round_proof.betas)
        comb_randomness = sponge.squeeze_f()
        new_folding_randomness = sponge.squeeze_f()

        scaling_factor = (
            verification_state.domain_size // self.parameters.folding_factor
        )
        num_repetitions = self.parameters.repetitions[verification_state.num_round]

        stir_randomness_indexes = sorted(
            list(
                set(sponge.squeeze_int(scaling_factor) for _ in range(num_repetitions))
            )
        )

        if not proof_of_work_verify(
            sponge,
            self.parameters.pow_bits[verification_state.num_round],
            round_proof.pow_nonce,
        ):
            return (False, verification_state)
        _ = sponge.squeeze_f()
        # Now, for each of the selected random points, we need to compute the folding of the
        # previous oracle

        oracle_answers = round_proof.queries_to_prev[0]
        folded_answers = self.compute_folded_evaluations(
            verification_state, stir_randomness_indexes, oracle_answers
        )
        return (True, verification_state)

    def compute_folded_evaluations(
        self,
        verification_state: VerificationState,
        stir_randomness_indexes: list,
        oracle_answers: list,
    ) -> list:
        return []
