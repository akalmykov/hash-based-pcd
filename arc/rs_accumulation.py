import math
from typing import List

import numpy as np
from stir.constants import FIELD192
import galois

from stir.domain import Domain
from stir.quotient import quotient_codeword
from stir.fiat_shamir import Blake3Sponge
from stir.merkle import MerkleTree
from arc.combine import combine_many, combine_poly


class PublicParameters:
    def __init__(self, GF: galois.GF, s: int, lambda_value: int):
        self.field = GF
        self.s = s  # out-of-domain repetitions
        self.t = None  # in-domain repetitions are calculated later
        self.lambda_value = lambda_value


class RSAccumulation:
    def __init__(self, p: PublicParameters, polys: List[galois.Poly]):
        self.p = p
        self.m = len(polys)
        self.polys = polys
        self.rate = 2
        self.d_max = max([len(poly) for poly in polys])
        # print("d_max = ", self.d_max)
        # print("rate = ", self.rate)
        self.L = Domain(self.p.field, self.d_max, self.rate)
        # print("L.size = ", self.L.size)
        assert self.L.size > self.d_max
        # print("|L| = ", self.L.size)
        self.rho = (self.d_max + 1) / self.L.size
        # print("(d_max+1) / |L| = ", self.rho)
        required_field_size = int(
            2**self.p.lambda_value * 10**7 * self.m * self.d_max**3 * self.rho ** (-3.5)
        )
        # print("required_field_size = ", required_field_size)
        # print("field size = ", self.p.field.characteristic)
        assert self.p.field.characteristic >= required_field_size
        self.eta = 0.01
        self.delta = 1 - math.sqrt(self.rho) - self.eta
        self.delta_bound = (
            1
            - 1.05 * math.sqrt(self.rho)
            - self.p.lambda_value / (self.L.size * (-math.log(1 - self.delta)))
        )
        # print("delta = ", self.delta)
        # print("delta_bound = ", self.delta_bound)
        min_distance = self.L.size - self.d_max + 1
        # print("minimum distance = ", min_distance)
        # print("relative minimum weight = ", min_distance / self.L.size)
        unique_decoding_distance = math.floor((min_distance - 1) / 2)
        # print("unique decoding distance = ", unique_decoding_distance)
        # print(
        #     "relative unique decoding distance = ",
        #     unique_decoding_distance / self.L.size,
        # )
        self.p.t = math.floor(self.p.lambda_value / (-math.log(1 - self.delta)))
        print("t = ", self.p.t)
        self.sponge = Blake3Sponge(self.p.field, 191)
        # TODO: absorb evals of polys into sponge

    def virtual_poly_quotient(self, r, points):
        poly = combine_poly(self.polys, r)
        evaluations = poly(points)
        ans_polynomial = galois.lagrange_poly(points, evaluations)
        x = galois.Poly.Identity(self.p.field)
        vanishing_polynomial = galois.Poly.One(self.p.field)
        for s in points:
            vanishing_polynomial *= x - s
        return (poly + ans_polynomial) // vanishing_polynomial

    def accumulate(self):
        elements = np.array(list(self.L.domain.elements()))
        r = self.sponge.squeeze_f()
        evals = combine_many(self.polys, r, elements)
        mt = MerkleTree([[e] for e in evals])
        self.sponge.absorb(bytes(mt.root()))
        x_out = self.sponge.squeeze_several_f(self.p.s)
        ood_evals = combine_many(self.polys, r, x_out)
        # TODO commit to ood_evals?
        x_in = []
        x_in_indexes = []
        while len(x_in) < self.p.t:
            x_in_index = self.sponge.squeeze_int(self.L.size)
            if x_in_index not in x_in_indexes:
                new_x_in = self.L.domain.element(x_in_index)
                assert new_x_in in elements, "x_in not in elements"
                x_in_indexes.append(x_in_index)
                x_in.append(new_x_in)
        poly_quotient = self.virtual_poly_quotient(r, self.p.field(x_out + x_in))
        fill = poly_quotient(x_in)

        # query phase
        x_in_indexes_sorted = sorted(x_in_indexes)
        f_oracle_evals = [[evals[i]] for i in x_in_indexes_sorted]
        proof = mt.generate_multi_proof(x_in_indexes_sorted)
        if not proof.verify(mt.root(), f_oracle_evals):
            print("mt verification failed")
            return
        else:
            print("mt verified")

        S = [x for x in x_out] + [x for x in x_in]
        ans_vals = [e for e in ood_evals] + [evals[i] for i in x_in_indexes]
        assert len(S) == len(ans_vals)
        S_to_ans = dict(zip([int(x) for x in S], ans_vals))

        self.d = self.d_max - self.L.size
        # c := Quotient(·, S, Ans, Fill) => c := Quotient(f, S, Ans, Fill)
        self.c = quotient_codeword(
            self.p.field, elements, evals, S, S_to_ans, poly_quotient
        )  # L


def test_rs_accumulation():
    GF = galois.GF(FIELD192)
    p = PublicParameters(GF, 1, 128)
    with open("./stir/test/original-stir-traces/poly_coeffs.json", "r") as f:
        import json

        fixed_poly = json.load(f)
    poly_coeffs = GF([int(fixed_poly[i]) for i in range(len(fixed_poly))])

    for _ in range(10):
        polys = []
        n = 100
        for _ in range(n):
            poly_coeffs = GF.Random(64, seed=1)
            polys.append(galois.Poly(poly_coeffs, field=GF, order="asc"))

        prover = RSAccumulation(p, polys)
        import time

        start_time = time.time()
        prover.accumulate()
        end_time = time.time()
        total_time = end_time - start_time
        print(f"rs accumulation time per MiMc instance: {total_time / n:.6f} seconds")
    prover = RSAccumulation(p, [galois.Poly(poly_coeffs, field=GF, order="asc")])
    prover.accumulate()
    # # Time the accumulate method over 2000 executions
    # import time

    # start_time = time.time()
    # for i in range(10):
    #     prover.accumulate()
    # end_time = time.time()

    # total_time = end_time - start_time
    # average_time = total_time / 10

    # print(f"Total time for 10 executions: {total_time:.4f} seconds")
    # print(f"Average time per accumulate(): {average_time:.6f} seconds")
    # print(f"Average time per accumulate(): {average_time * 1000:.3f} milliseconds")
