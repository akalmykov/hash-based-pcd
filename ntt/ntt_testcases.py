from stir.domain import Domain
from stir.constants import FIELD192
import galois
import json


def test_ntt():
    GF = galois.GF(FIELD192)
    domain = Domain(GF, 2**10, 2)
    evals = []
    with open("./stir/test/original-stir-traces/evals.json", "r") as f:
        evals = json.load(f)
    t_evals = GF([int(evals[i]) for i in range(len(evals))])

    poly_coeffs = []
    with open("./stir/test/original-stir-traces/poly_coeffs.json", "r") as f:
        poly_coeffs = json.load(f)
    t_poly_coeffs = GF([int(poly_coeffs[i]) for i in range(len(poly_coeffs))])

    my_evals = domain.evaluate_poly_coeff(t_poly_coeffs)
    assert (my_evals == t_evals).all()

    test_case = {
        "degree": 2**10,
        "log_rho_inv": 2,
        "poly_coeffs": t_poly_coeffs,
        "evals": t_evals,
    }
    return test_case


if __name__ == "__main__":
    ntt_test_cases = test_ntt()
