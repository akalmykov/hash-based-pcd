import galois

from stir.constants import FIELD192
from stir.domain import Domain
from dataclasses import dataclass
from typing import List


@dataclass
class Racc:
    e: int  # field element (scalar)
    v: List[galois.Array]  # vector of field elements
    f: List[galois.Array]  # vector of field elements


def test_r1cs_to_racc():
    import json

    GF = galois.GF(FIELD192)

    # read r1cs from json
    with open("./mimc_r1cs.json", "r") as f:
        r1cs = json.load(f)

    w = []
    for w_str in r1cs["witness"]:
        w.append(GF(int(w_str, 16)))
    m = len(w)
    # print("m", m)
    domain = Domain(GF, len(w), 2)
    # print("domain.size", domain.size)
    # print("domain.get_size()", domain.get_size())
    w_encoded = domain.evaluate_poly_coeff(w)
    # print("len(w_encoded)", len(w_encoded))
    racc = Racc(0, [GF.Random() for _ in range(1, m)], w_encoded)
