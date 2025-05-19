from stir.prover import Prover
import json

if __name__ == "__main__":
    prover = Prover()
    with open("./stir/test/original-stir-traces/poly_coeffs.json", "r") as f:
        fixed_poly = json.load(f)
    poly_coeffs = prover.field([int(fixed_poly[i]) for i in range(len(fixed_poly))])

    prover.commit(poly_coeffs)
    prover.prove()
