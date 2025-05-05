import galois


def stack_evals(evals, folding_factor):
    assert len(evals) % folding_factor == 0
    size_of_the_new_domain = len(evals) // folding_factor
    stacked_evals = []
    for i in range(size_of_the_new_domain):
        new_evals = []
        for j in range(folding_factor):
            new_evals.append(evals[i + size_of_the_new_domain * j])
        stacked_evals.append(new_evals)
    return stacked_evals


def is_power_of_two(n):
    return n & (n - 1) == 0


def next_power_of_two(n):
    return 1 << (n - 1).bit_length()


def get_two_adicity(GF):
    n = GF.order - 1

    # Count trailing zeros to get two-adicity
    adicity = 0
    while n & 1 == 0:
        adicity += 1
        n >>= 1

    return adicity
