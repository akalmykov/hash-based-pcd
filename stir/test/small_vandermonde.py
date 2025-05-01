import numpy as np


def find_16th_root_of_unity_in_F17():
    """Find a primitive 16th root of unity in the finite field F17."""
    # We need to find an element ω such that ω^16 ≡ 1 (mod 17)
    # and ω^k ≠ 1 (mod 17) for 1 ≤ k < 16

    for candidate in range(1, 17):
        # Check if candidate^16 ≡ 1 (mod 17)
        if pow(candidate, 16, 17) == 1:
            # Check if it's primitive (i.e., not a root of unity of lower order)
            is_primitive = True
            for k in range(1, 16):
                if pow(candidate, k, 17) == 1:
                    is_primitive = False
                    break

            if is_primitive:
                return candidate

    return None


# Let's try another approach. Since the multiplicative group of F17 is of order 16,
# we can find a generator of the group. Any generator raised to the appropriate power
# can give us a primitive 16th root of unity.


def find_generator_of_F17():
    """Find a generator of the multiplicative group of F17."""
    # The multiplicative group of F17 has order 16
    # A generator g must satisfy that the smallest n where g^n ≡ 1 (mod 17) is 16

    for g in range(2, 17):
        # Calculate orders of elements
        orders = set()
        for i in range(1, 17):
            orders.add(pow(g, i, 17))

        # If we have 16 distinct elements, then g is a generator
        if len(orders) == 16:
            return g

    return None


def create_vandermonde_matrix(omega, size, mod):
    """Create a Vandermonde matrix using omega in the finite field F_mod."""
    matrix = []
    for i in range(size):
        row = []
        for j in range(size):
            # Calculate omega^(i*j) mod 17
            value = pow(omega, i * j, mod)
            row.append(value)
        matrix.append(row)

    return matrix


# Find a generator of F17
generator = find_generator_of_F17()
print(f"Generator of F17: {generator}")

# To get a primitive 16th root of unity, we need generator^1
omega = generator
print(f"A primitive 16th root of unity in F17: {omega}")

# Verify
powers = [pow(omega, i, 17) for i in range(16)]
print(f"Powers of omega: {powers}")
print(f"Check omega^16 ≡ 1 (mod 17): {pow(omega, 16, 17) == 1}")

# Create a 4×4 Vandermonde matrix using this root of unity
size = 4
vandermonde_matrix = create_vandermonde_matrix(omega, size, 17)

print(f"\n{size}×{size} Vandermonde matrix using the 16th root of unity in F17:")
for row in vandermonde_matrix:
    print(" ".join(f"{val:2}" for val in row))

# Create the full 16×16 Vandermonde matrix
full_size = 16
full_vandermonde_matrix = create_vandermonde_matrix(omega, full_size, 17)

print(f"\nFull {full_size}×{full_size} Vandermonde matrix:")
for row in full_vandermonde_matrix:
    print(" ".join(f"{val:2}" for val in row))
