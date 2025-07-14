import galois

# Define base field F_q
q = 2**61 - 1
Fq = galois.GF(q)

# Find irreducible polynomial of degree 3 over F_q
f = galois.irreducible_poly(q, degree=2)
print("Irreducible polynomial of degree 2 over F_q:", f)
# Define the field F_q[x]/(f)
Fq2 = galois.GF(q**2, irreducible_poly=f)

log_n = 4
n = 2**log_n
v = Fq.Random(n)
galois.Poly
