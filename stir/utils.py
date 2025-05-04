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
