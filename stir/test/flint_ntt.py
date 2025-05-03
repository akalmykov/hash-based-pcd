from flint import fmpz_mod_poly_ctx, fmpz_mod_poly
import random

# Use 998244353 — a well-known NTT-friendly prime
p = 998244353
ctx = fmpz_mod_poly_ctx(p)

N = 2**20

# Create polynomials using proper FLINT syntax
a = fmpz_mod_poly([1, 2, 3], ctx)  # Initialize with coefficients list
b = fmpz_mod_poly([4, 5, 6], ctx)  # Initialize with coefficients list

# Multiply using FLINT — internally uses NTT
c = a * b

# Print results
print("Polynomial a:", a)
print("Polynomial b:", b)
print("Product c = a * b:", c)
