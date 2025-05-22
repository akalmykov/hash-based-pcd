from .constants import FIELD192
from .prover import Prover, Proof
from .verifier import Verifier
import galois
from .merkle import MerkleTree
from .stir_parameters import Parameters, DEFAULT_PARAMETERS


def run_stir(prove: bool, verify: bool):
    gf = galois.GF(FIELD192)
    if prove:
        prover = Prover(gf)
        with open("./stir/test/original-stir-traces/poly_coeffs.json", "r") as f:
            import json

            fixed_poly = json.load(f)
        poly_coeffs = prover.field([int(fixed_poly[i]) for i in range(len(fixed_poly))])

        prover.commit(poly_coeffs)
        prover.merkle_tree.serialize("./stir/test/merkle_tree.pkl")

        proof = prover.prove()
        proof.serialize("./stir/test/proof.pkl")
    else:
        print("Skipped proving")
    if verify:

        mt = MerkleTree.load("./stir/test/merkle_tree.pkl")
        proof = Proof.load("./stir/test/proof.pkl")

        verifier = Verifier(gf, proof.parameters)
        assert verifier.verify(mt, proof)
    else:
        print("Skipped verification")
