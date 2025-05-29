from .prover import Proof
from .stir_parameters import FullParameters
from .merkle import MerkleTree
from .fiat_shamir import Blake3Sponge
from .domain import Domain, EvaluationDomainConfig, Radix2EvaluationDomain
from dataclasses import dataclass
from typing import Any
import galois
from enum import Enum
from typing import Tuple
from .prover import RoundProof
from .pow import proof_of_work_verify
from .interpolation import fft_interpolate
import numpy as np


class BaseOracle:
    def __init__(self, GF):
        self.GF = GF

    def common_factor_scale(self) -> galois.FieldArray:
        raise NotImplementedError("This method should be overridden in subclasses.")

    def get_denominators(self, query_set: list) -> list:
        raise NotImplementedError("This method should be overridden in subclasses.")


class InitialOracle(BaseOracle):
    def __init__(self, GF):
        self.GF = GF

    def common_factor_scale(self) -> galois.FieldArray:
        return self.GF(0)

    def get_denominators(self, query_set: list) -> list:
        return self.GF([1] * len(query_set))


class VirtualOracle(InitialOracle):
    def __init__(self, GF, comb_randomness, interpolating_polynomial, quotient_set):
        super().__init__(GF)
        self.comb_randomness = comb_randomness
        self.interpolating_polynomial = interpolating_polynomial
        self.quotient_set = quotient_set

    def common_factor_scale(self):
        return self.comb_randomness

    def get_denominators(self, query_set: list) -> list:
        return [
            np.prod([eval_point - x for x in self.quotient_set])
            for eval_point in query_set
        ]


@dataclass
class VerificationState:
    oracle: BaseOracle
    domain_gen: galois.FieldArray
    domain_size: int
    domain_offset: galois.FieldArray
    root_of_unity: galois.FieldArray
    folding_randomness: galois.FieldArray
    num_round: int

    def query(
        self,
        evaluation_point: galois.FieldArray,
        value_of_prev_oracle: galois.FieldArray,
        common_factors_inverse: galois.FieldArray,
        denom_hint: galois.FieldArray,
        ans_eval: galois.FieldArray,
    ) -> galois.FieldArray:
        if isinstance(self.oracle, InitialOracle):
            return value_of_prev_oracle
        elif isinstance(self.oracle, VirtualOracle):
            num_terms = len(self.oracle.quotient_set)

            # Compute quotient evaluation
            quotient_evaluation = (value_of_prev_oracle - ans_eval) * denom_hint

            common_factor = evaluation_point * self.oracle.comb_randomness

            if common_factor != self.oracle.GF(1):
                scale_factor = (
                    self.oracle.GF(1) - common_factor ** (num_terms + 1)
                ) * common_factors_inverse
            else:
                scale_factor = self.oracle.GF(num_terms + 1)

            return quotient_evaluation * scale_factor


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
            oracle=InitialOracle(self.GF),
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
        scaling_factor = (
            verification_state.domain_size // self.parameters.folding_factor
        )
        generator = verification_state.domain_gen**scaling_factor
        coset_offsets = [
            verification_state.domain_offset * (verification_state.domain_gen**i)
            for i in stir_randomness_indexes
        ]
        scales: list[galois.FieldArray] = []
        scale: galois.FieldArray = self.GF(1)
        for _ in range(self.parameters.folding_factor):
            scales.append(scale.copy())  # otherwise, it will be a reference
            scale *= generator

        query_sets = [
            [coset_offset * scales[j] for j in range(self.parameters.folding_factor)]
            for coset_offset in coset_offsets
        ]

        common_factor_scale = verification_state.oracle.common_factor_scale()
        global_common_factors = [
            [self.GF(1) - common_factor_scale * x for x in query_set]
            for query_set in query_sets
        ]

        global_denominators = [
            verification_state.oracle.get_denominators(query_set)
            for query_set in query_sets
        ]
        size = self.parameters.folding_factor
        to_invert = []
        global_common_factors_len = len(global_common_factors)
        for common_factors in global_common_factors:
            to_invert.extend(common_factors)
        for denominators in global_denominators:
            to_invert.extend(denominators)
        to_invert.extend(coset_offsets)
        to_invert.append(generator)
        to_invert.append(size)
        to_invert = self.GF(to_invert) ** -1

        # Extract the inverse values
        size_inv = to_invert[-1]
        generator_inv = to_invert[-2]
        coset_offsets_inv = to_invert[-len(coset_offsets) - 2 : -2]

        # Extract remaining elements
        remaining = to_invert[: -len(coset_offsets) - 2]

        # Chunk into lists of size folding_factor
        chunked = [
            remaining[i : i + self.parameters.folding_factor]
            for i in range(0, len(remaining), self.parameters.folding_factor)
        ]

        # Split into common_factors_inv and denominators_inv
        common_factors_inv = chunked[:global_common_factors_len]
        denominators_inv = chunked[global_common_factors_len:]

        # Compute evaluations of answers for each coset
        evaluations_of_ans = []
        for coset_offset, coset_offset_inv in zip(coset_offsets, coset_offsets_inv):
            if isinstance(verification_state.oracle, InitialOracle):
                evaluations_of_ans.append(self.GF([1] * self.parameters.folding_factor))
            elif isinstance(verification_state.oracle, VirtualOracle):
                offset_pow_size = coset_offset**self.parameters.folding_factor

                # Set up parameters for polynomial evaluation
                domain_params = EvaluationDomainConfig(
                    size=self.parameters.folding_factor,
                    size_as_field_element=size,
                    size_inv=size_inv,
                    group_gen=generator,
                    group_gen_inv=generator_inv,
                    offset=coset_offset,
                    offset_inv=coset_offset_inv,
                    offset_pow_size=offset_pow_size,
                )
                # Evaluate the interpolating polynomial over the domain
                virtual_function = verification_state.oracle
                evals = Radix2EvaluationDomain.from_config(
                    self.GF, domain_params
                ).evaluate_poly_coeff(virtual_function.interpolating_polynomial)
                evaluations_of_ans.append(evals)
        scaled_offset = verification_state.domain_offset**self.parameters.folding_factor
        folded_answers = []
        print("stir_randomness_indexes length:", len(stir_randomness_indexes))
        print("coset_offsets length:", len(coset_offsets))
        print("coset_offsets_inv length:", len(coset_offsets_inv))
        print("query_sets length:", len(query_sets))
        print("common_factors_inv length:", len(common_factors_inv))
        print("denominators_inv length:", len(denominators_inv))
        print("evaluations_of_ans length:", len(evaluations_of_ans))

        for i, (
            stir_randomness_index,
            coset_offset,
            coset_offset_inv,
            query_set,
            common_factors_inv,
            denominators_inv,
            evaluation_of_ans,
        ) in enumerate(
            zip(
                stir_randomness_indexes,
                coset_offsets,
                coset_offsets_inv,
                query_sets,
                common_factors_inv,
                denominators_inv,
                evaluations_of_ans,
            )
        ):
            stir_randomness = scaled_offset * verification_state.domain_gen ** (
                self.parameters.folding_factor * stir_randomness_index
            )
            f_answers = [
                verification_state.query(
                    x,
                    oracle_answers[i][j],
                    common_factors_inv[j],
                    denominators_inv[j],
                    evaluation_of_ans[j],
                )
                for j, x in enumerate(query_set)
            ]

            folded_answer_poly = fft_interpolate(
                self.GF,
                generator,
                coset_offset,
                generator_inv,
                coset_offset_inv,
                size_inv,
                f_answers,
            )
            if i == 0:
                print("generator:", generator)
                print("coset_offset:", coset_offset)
                print("generator_inv:", generator_inv)
                print("coset_offset_inv:", coset_offset_inv)
                print("size_inv:", size_inv)
                print("\n***** folded_answer_poly", folded_answer_poly)
            # for x in f_answers:
            #     print(int(x), end=",")

            folded_answer_eval = folded_answer_poly(
                verification_state.folding_randomness
            )
            folded_answers.append((stir_randomness, folded_answer_eval))

        return []
