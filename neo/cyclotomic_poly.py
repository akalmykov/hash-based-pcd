from sage.all import *


# m is the cyclotomic polynomial index
def tau(m):
    return m if (m % 2) != 0 else m / 2


# [LS18, eprint 2017-523] Thm 1.1, pg 4
# m is the cyclotomic polynomial index
# p is the prime
# z is any divisor of m
# This tests for the condition for thm 1.1 to hold
def thm1_1_cond(m, p, z):
    cond1 = (p % z) == 1
    cond2 = Mod(p, m).multiplicative_order() == m / z
    return cond1 and cond2


# [LS18, eprint 2017-523] Thm 1.1, pg 4
# p is the prime
# z is any divisor of m
# lInf bound for elements to be invertible
# given that m,p,z satisfy thm 1.1 cond
def thm1_1_inv_bound(p, z):
    return (1 / s1(z) * p ^ (1 / euler_phi(z))).n()


def thm1_1_num_factors(z):
    return euler_phi(z)


# Output divisors of m
def divisors(m):
    zs = list()
    for i in range(1, m + 1):
        if m % i == 0:
            zs.append(i)
    return zs


# [LS18, eprint 2017-523] pg 6, pg 9
# We only consider prime power cyclotomics
# m is the cyclotomic polynomial index
def s1(m):
    return sqrt(tau(m))


# checks if cyclotomic index m is power of two
def is_pow2(m):
    return sum(m.digits(2)) == 1


# [MR09] lattice-based cryptography


# makes sure characteristic does not lead
# to trivial bound
def non_trivial(q, n, d, delta):
    return (q / 2).n() >= (2 ^ (2 * sqrt(n * d * log(q, 2) * log(delta, 2)))).n()


# [AL21] eprint Prop 2. 2021/202
# for all u,v in R, |u*v| / |v| <= gamma*|u|
# outputs T = gamma * |u|
# assumes we are only testing prime powers
def expansion_factor(m, norm):
    if is_pow2(m):
        return euler_phi(m) * norm
    else:
        return 2 * euler_phi(m) * norm


# p is prime
# max_idx is max cyclotomic index
# outputs list of (m, z)
def candidates(p, min_idx=10, max_idx=200):
    # prime powers
    possible_indices = [i for i in range(min_idx, max_idx) if len(factor(i)) == 1]
    c = list()
    for m in possible_indices:
        zs = divisors(m)
        for z in zs:
            if thm1_1_cond(m, p, z):
                c.append((Integer(m), Integer(z)))
    return c


def pre_filter(q, cyclotomic_index, z, n, m, chals):
    chals_norm = max({abs(c) for c in chals})
    chals_max_diff = chals[-1] - chals[0]
    delta = 1.0045  # root hermite factor, chosen from [ESSLL19] eprint 2018/773
    phi = cyclotomic_polynomial(cyclotomic_index)  # index cyclotomic polynomial
    d = phi.degree()  # degree of cyclotomic
    # return non_trivial(q, n, d, delta) and chals_max_diff < thm1_1_inv_bound(q, z) and log(len(chals)^d,2).n() >= 120
    # We remove non_trivial(...) because we use the lattice estimator for hardness
    return chals_max_diff < thm1_1_inv_bound(q, z) and log(len(chals) ^ d, 2).n() >= 120


def info(q, cyclotomic_index, z, n, m, chals):
    chals_norm = max({abs(c) for c in chals})
    chals_max_diff = chals[-1] - chals[0]
    delta = 1.0045  # root hermite factor, chosen from [ESSLL19] eprint 2018/773
    phi = cyclotomic_polynomial(cyclotomic_index)  # index cyclotomic polynomial
    d = phi.degree()  # degree of cyclotomic
    T = expansion_factor(cyclotomic_index, chals_norm)
    # Bounds for MSIS to be hard
    # [MR09] lattice-based cryptography pg 6
    # [CMNW24] pg 38 eprint 2024/281
    MSIS_B_l2_bound = min(q, 2 ^ (2 * sqrt(n * d * log(q, 2) * log(delta, 2))))
    MSIS_B_linf_bound = MSIS_B_l2_bound / sqrt(m * d)
    # We need MSIS infinity bound 8TB to be hard
    B = MSIS_B_linf_bound / (8 * T)
    print("####")
    print("Cyclotomic idx:", cyclotomic_index)
    print("Cyclotomic Poly:", phi)
    print("z:", z)
    # print("Prime is non-trivial?", non_trivial(q, n, d, delta))
    print("Csmall norm is small enough?", chals_max_diff < thm1_1_inv_bound(q, z))
    print("Csmall large enough?", log(len(chals) ^ d, 2).n() >= 120)
    print("Degree of Cyclotomic:", d)
    # print("log(B):", log(B, 2).n())
    print("Expansion Factor T:", T)
    print("Invertible Norm bound:", thm1_1_inv_bound(q, z))
    print("log(|C_Small|):", log(len(chals) ^ d, 2).n())
    print("Factors of Cyclotomic:", thm1_1_num_factors(z))
    print()


def possible_settings(q, n, m, chals):
    for cyclotomic_index, z in candidates(q):
        if pre_filter(q, cyclotomic_index, z, n, m, chals):
            info(q, cyclotomic_index, z, n, m, chals)
        else:
            delta = 1.0045
            d = cyclotomic_polynomial(cyclotomic_index).degree()
            print(
                "[Does not satisfy security requirements] index: {}, degree: {}, z: {}, non_trivial"
            )


if __name__ == "__main__":
    # Primes:
    GL = 2 ^ 64 - 2 ^ 32 + 1
    AGL = GL - 32
    print("###############################")
    print("AGL ###############################")
    print("###############################")
    # MSIS settings
    n = 13  # rows, kappa in latticefold
    m = 2 ^ 26  # cols
    # Small Challenge set
    chals = [-1, 0, 1, 2]
    possible_settings(AGL, n, m, chals)
    print("###############################")
    print("M61 ###############################")
    print("###############################")
    # MSIS settings
    n = 16  # rows, kappa in latticefold
    m = 2 ^ 22  # cols
    # Small Challenge set
    chals = [-2, -1, 0, 1, 2]
    possible_settings(2 ^ 61 - 1, n, m, chals)
    print("###############################")
    print("GL ###############################")
    print("###############################")
    # MSIS settings
    n = 16  # rows, kappa in latticefold
    m = 2 ^ 24  # cols
    # Small Challenge set
    chals = [-2, -1, 0, 1, 2]
    possible_settings(GL, n, m, chals)
    print("###############################")
