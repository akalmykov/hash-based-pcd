import galois

from constants import TEST_FIELD


def poly_quotient(poly, points):
    GF = type(points[0])
    evaluations = poly(points)
    ans_polynomial = galois.lagrange_poly(points, evaluations)
    x = galois.Poly.Identity(GF)
    vanishing_polynomial = galois.Poly.One(GF)
    for s in points:
        vanishing_polynomial *= x - s
    return (poly + ans_polynomial) // vanishing_polynomial


# This is almost Quotient(f, S, Ans, Fill) in the paper
# If x ∈ S, the value is given by Fill(x)
# If x ∉ S, the value is given by (f(x) - AnsPoly(x)) / VanishingPoly(x)
# where VanishingPoly(x) = 0 for x ∈ S
def quotient(claimed_eval, evaluation_point, answers: list[tuple[int, int]]):
    GF = type(claimed_eval)
    if evaluation_point in [answer[0] for answer in answers]:
        # return answers[evaluation_point]
        raise ValueError("Evaluation point is in the domain")
    else:
        ans_polynomial = galois.lagrange_poly(
            GF([answer[0] for answer in answers]),
            GF([answer[1] for answer in answers]),
        )
        ans_eval = ans_polynomial(evaluation_point)
        numerator = claimed_eval - ans_eval
        denominator = 1
        for answer in answers:
            denominator *= evaluation_point - answer[0]
        return numerator / denominator


if __name__ == "__main__":
    #   ---- poly_utils::quotient::tests::test_quotient stdout ----
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
    # 11171883536905842650 * x^10
    # points: [0, 1]
    # quotient_poly:
    # 14814265841813440540 +
    # 8761597112993598527 * x +
    # 11723079719027242991 * x^2 +
    # 2065466484697509892 * x^3 +
    # 8751868778253374432 * x^4 +
    # 963294955979777575 * x^5 +
    # 623279022173724428 * x^6 +
    # 14334562976994545075 * x^7 +
    # 11171883536905842650 * x^8
    # ans: [(0, 11610679913765557349), (1, 15116961601664466783)]
    # test_point: 6160416209525857734
    # quotient_poly.evaluate(&test_point): 5912928968069917069

    gf64 = galois.GF(TEST_FIELD)
    degree = 10
    poly = galois.Poly(gf64.Random(degree + 1), field=gf64)
    # poly = galois.Poly(
    #     list(
    #         reversed(
    #             [
    #                 11610679913765557349,
    #                 7138759915500053215,
    #                 6052668728819842013,
    #                 15485261463380939857,
    #                 9657613234329733099,
    #                 11760341775858719781,
    #                 7788573822273596857,
    #                 340015933806053147,
    #                 4735460114593763674,
    #                 3162679440088702425,
    #                 11171883536905842650,
    #             ]
    #         )
    #     ),
    #     field=gf64,
    # )
    points = gf64([0, 1])
    quotient_poly = poly_quotient(poly, points)
    # print(quotient)
    ans = list(zip(points, poly(points)))
    # print(ans)
    # test_point = gf64(6160416209525857734)
    test_point = gf64.Random(1)[0]
    # print(quotient_poly(test_point))
    # print(quotient(poly(test_point), test_point, ans))
    assert quotient(poly(test_point), test_point, ans) == quotient_poly(test_point)
