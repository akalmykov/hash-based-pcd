import numpy as np

from stir.domain import Radix2EvaluationDomain, EvaluationDomainConfig, Domain
import galois

from stir.utils import stack_evals
from stir.utils import is_power_of_two


def fft_interpolate_simple(
    GF, generator, coset_offset, points
):
    points = list(points)
    folding_factor = len(points)
    assert is_power_of_two(folding_factor), "Folding factor must be a power of two"
    size_as_field_element = GF(folding_factor)
    size_inv = size_as_field_element ** (-1)
    coset_offset_inv = coset_offset ** (-1)
    print("coset_offset_inv", coset_offset_inv)
    generator_inv = generator ** (-1)

    config = EvaluationDomainConfig(
        size=folding_factor,
        size_as_field_element=size_as_field_element,
        size_inv=size_inv,
        group_gen=generator,
        group_gen_inv=generator_inv,
        offset=coset_offset,
        offset_inv=coset_offset_inv,
        offset_pow_size=coset_offset**folding_factor,
    )

    # Create domain from config
    domain = Radix2EvaluationDomain.from_config(GF, config)

    # Perform inverse NTT to get polynomial coefficients
    if domain.offset != GF(1):
        alpha_inv_powers = domain.offset_inv ** (np.arange(folding_factor))
        coeffs = galois.intt(points, domain.size, GF.order) * alpha_inv_powers
    else:
        coeffs = galois.intt(points, domain.size, GF.order)
    return galois.Poly(list(reversed(coeffs)), field=GF)

def fft_interpolate(
    GF, generator, coset_offset, generator_inv, coset_offset_inv, size_inv, points
):
    # Convert points to list if needed
    points = list(points)
    folding_factor = len(points)

    assert is_power_of_two(folding_factor), "Folding factor must be a power of two"

    size_as_field_element = GF(folding_factor)

    # Create configuration for the domain
    config = EvaluationDomainConfig(
        size=folding_factor,
        size_as_field_element=size_as_field_element,
        size_inv=size_inv,
        group_gen=generator,
        group_gen_inv=generator_inv,
        offset=coset_offset,
        offset_inv=coset_offset_inv,
        offset_pow_size=coset_offset**folding_factor,
    )

    # Create domain from config
    domain = Radix2EvaluationDomain.from_config(GF, config)

    # Perform inverse NTT to get polynomial coefficients
    if domain.offset != GF(1):
        alpha_inv_powers = domain.offset_inv ** (np.arange(folding_factor))
        coeffs = galois.intt(points, domain.size, GF.order) * alpha_inv_powers
    else:
        coeffs = galois.intt(points, domain.size, GF.order)

    # Return polynomial with coefficients in standard order (highest degree first)
    return galois.Poly(list(reversed(coeffs)), field=GF)


if __name__ == '__main__':
    from constants import TEST_FIELD

    GF = galois.GF(TEST_FIELD)
    degree = 15
    folding_factor = 8
    polynomial = galois.Poly([1]*(degree+1), field=GF)
    print(polynomial)
    domain = Domain(GF, degree, 3)
    evaluations = domain.evaluate_poly(polynomial)

    elements = list(domain.domain.elements())

    reshaped_elements = stack_evals(elements, folding_factor)
    reshaped_evaluations = stack_evals(evaluations, folding_factor)

    def print_gf_array(arr):
        for elem in arr:
            print(elem, end=",")
        print()

    g_evaluations = []
    for i, evals in enumerate(reshaped_evaluations):
        x = []
        y = []
        for j, val in enumerate(evals):
            x.append(reshaped_elements[i][j])
            y.append(val)
        g_evaluations.append(galois.lagrange_poly(GF(x), GF(y)))

    print("g_evaluations", g_evaluations)
    generator = domain.element(domain.domain.size // folding_factor)
    print("generator", generator)

    fft_polys = []
    for i, evals in enumerate(reshaped_evaluations):
        coset_offset = domain.domain.element(i)
        print("coset_offset", coset_offset)
        fft_polys.append(fft_interpolate_simple(GF, generator, coset_offset, evals))

    assert fft_polys == g_evaluations

    generator_inv = generator ** (-1)
    size_as_field_element = GF(folding_factor)
    size_inv = size_as_field_element ** (-1)


    fft_polys2 = []
    for i, evals in enumerate(reshaped_evaluations):
        coset_offset = domain.domain.element(i)
        coset_offset_inv = coset_offset ** (-1)
        print("coset_offset", coset_offset)
        fft_polys2.append(fft_interpolate(GF, generator, coset_offset, generator_inv, coset_offset_inv, size_inv, evals))
    assert fft_polys2 == g_evaluations
