import numpy as np
import galois


class BivariatePolynomial:
    def __init__(self, matrix: galois.FieldArray):
        self.matrix = matrix

    def degree_x(self):
        return self.rows() - 1

    def rows(self):
        return len(self.matrix)

    def degree_y(self):
        return self.cols() - 1

    def cols(self):
        return len(self.matrix[0])

    def evaluate(self, x, y):
        field = type(x)
        result = field(0)  # Create a zero element in the same Galois Field
        for row in range(self.rows()):
            for col in range(self.cols()):
                term = self.matrix[row][col] * x**row * y**col
                result = result + term
        return result

    def fold_by_col(self, alpha):
        field = type(alpha)
        transposed = self.matrix.T

        result = field([0 for _ in range(self.rows())])

        pow_alpha = 1
        for col in transposed:
            scaled_col = col * pow_alpha
            result = result + scaled_col
            pow_alpha = pow_alpha * alpha

        return galois.Poly(list(reversed(result)), field=field)


def to_coefficient_matrix(poly_coeffs, rows, cols):
    if len(poly_coeffs) > rows * cols:
        raise ValueError("Degree of polynomial is too large for matrix")

    gf = type(poly_coeffs[0])
    matrix = gf([[0 for _ in range(cols)] for _ in range(rows)])

    for i, coeff in enumerate(poly_coeffs):
        row = i // cols
        col = i % cols
        matrix[row][col] = coeff

    return BivariatePolynomial(matrix)


if __name__ == "__main__":
    gf64 = galois.GF(18446744069414584321)
    # The order of the coefficients is X^0, X^1, X^2, ..., X^d
    coeff = [0, 1, 2, 3, 4, 5]  # 0X^0 + 1X^1 + 2X^2 + 3X^3 + 4X^4 + 5X^5
    poly = gf64(coeff)
    # galois.Poly EXPECTS IN THE ORDER X^d, X^d-1, ..., X^0
    # coefficients ${a_{d}, a_{d-1}, \dots, a_1, a_0}$
    gpoly = galois.Poly(list(reversed(coeff)), field=gf64)
    print(gpoly.degree)
    rows = 3
    cols = 2
    bivariate_poly = to_coefficient_matrix(poly, rows, cols)
    print(bivariate_poly.matrix)
    assert bivariate_poly.matrix.shape == (rows, cols)
    for i in range(bivariate_poly.rows()):
        for j in range(bivariate_poly.cols()):
            assert bivariate_poly.matrix[i, j] == poly[i * bivariate_poly.cols() + j]
    random_points = gf64.Random(10)
    for point in random_points:
        eval1 = bivariate_poly.evaluate(point**cols, point)
        eval2 = gpoly(point)
        assert eval1 == eval2
    # 0 1  ==> 0x^0 * y^0 + 1x^0 * y^1  ===  0x^0 + 1x^1
    # 2 3  ==> 2x^1 * y^0 + 3x^1 * y^1  ===  2x^2 + 3x^3  ==> fold these columns by alpha
    # 4 5  ==> 4x^2 * y^0 + 5x^2 * y^1  ===  4x^4 + 5x^5
    # The relationship between the original polynomial and the bivariate one is:
    # f(z) = q(z^n, z)
    # y=x^cols=x^2
    # 0X^0 + 1X^1 + 2X^2 + 3X^3 + 4X^4 + 5X^5
    # https://mhostetter.github.io/galois/latest/basic-usage/poly-arithmetic/#special-arithmetic
