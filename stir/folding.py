import math
from bivariate import to_coefficient_matrix
import galois


def poly_fold(poly_coeffs, folding_factor, folding_randomness):

    q_poly = to_coefficient_matrix(
        poly_coeffs, math.ceil(len(poly_coeffs) / folding_factor), folding_factor
    )
    q_poly = q_poly.fold_by_col(folding_randomness)
    return q_poly


if __name__ == "__main__":

    #   ---- poly_utils::folding::tests::test_folding stdout ----
    # poly:
    # 11610679913765557349 +
    # 7138759915500053215 * x +
    # 6052668728819842013 * x^2 +
    # 15485261463380939857 * x^3 +
    # 9657613234329733099 * x^4 +
    # 11760341775858719781 * x^5 +
    # 7788573822273596857 * x^6 +
    # 340015933806053147 * x^7 +
    # 4735460114593763674 * x^8 +
    # 3162679440088702425 * x^9 +
    # 11171883536905842650 * x^10 +
    # 6160416209525857734 * x^11 +
    # 15418064586095432452 * x^12 +
    # 18357664592583857644 * x^13 +
    # 2365385105606375019 * x^14 +
    # 16411461134942939508 * x^15 +
    # 996908178701015277 * x^16
    # poly_fold:
    # 10410991352436654782 +
    # 9691999768066204014 * x +
    # 13119089905379579041 * x^2 +
    # 9488653491303862592 * x^3 +
    # 2102113245622691478 * x^4 +
    # 5080476445705962678 * x^5 +
    # 14972667201941799067 * x^6 +
    # 10635714502662735275 * x^7 +
    # 996908178701015277 * x^8
    # root_of_unity: 13797081185216407910
    # evalpoint: 17870292113338400769
    # beta_l: [13797081185216407910, 4649662884198176411]
    # f_answers: [(13797081185216407910, 2035048425270305004), (4649662884198176411, 18027330255546616611)]
    # folded: 16849275407437574520
    # poly_fold.evaluate(&evalpoint): 16849275407437574520
    gf64 = galois.GF(18446744069414584321)
    degree = 16
    # poly = gf64(
    # [
    # 11610679913765557349,
    # 7138759915500053215,
    # 6052668728819842013,
    # 15485261463380939857,
    # 9657613234329733099,
    # 11760341775858719781,
    # 7788573822273596857,
    # 340015933806053147,
    # 4735460114593763674,
    # 3162679440088702425,
    # 11171883536905842650,
    # 6160416209525857734,
    # 15418064586095432452,
    # 18357664592583857644,
    # 2365385105606375019,
    # 16411461134942939508,
    # 996908178701015277,
    # ]
    # )
    poly = gf64.Random(degree + 1, seed=1)
    folding_factor = 2
    folding_randomness = gf64(5)
    folded_poly = poly_fold(poly, folding_factor, folding_randomness)
    folded_poly = galois.Poly(list(reversed(folded_poly)), field=gf64)
    print(folded_poly)
    root_of_unity = gf64.primitive_root_of_unity(256)
    # print(root_of_unity)
    eval_point = root_of_unity**folding_factor
    # print(eval_point)
    beta_0 = root_of_unity
    beta_1 = root_of_unity ** (1 + 128)
    # print(beta_0)
    # print(beta_1)
    assert eval_point == beta_0**2 and eval_point == beta_1**2
    gpoly = galois.Poly(list(reversed(poly)), field=gf64)
    lagrange_poly = galois.lagrange_poly(
        gf64([beta_0, beta_1]), gf64([gpoly(beta_0), gpoly(beta_1)])
    )
    print("folded:", folded_poly(eval_point))
    assert lagrange_poly(folding_randomness) == folded_poly(eval_point)
    # https://aszepieniec.github.io/stark-anatomy/fri
