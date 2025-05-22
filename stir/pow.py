def proof_of_work(sponge, proof_of_work_bits: int) -> int | None:
    assert proof_of_work_bits <= 32, "proof_of_work_bits must be <= 32"

    if proof_of_work_bits == 0:
        return None

    nonce = 0
    while True:
        new_sponge = sponge.clone()
        nonce_bytes = nonce.to_bytes(length=8, byteorder="little")
        new_sponge.absorb(nonce_bytes)
        pow_bytes = new_sponge.squeeze_bytes(4)
        pow = int.from_bytes(pow_bytes, byteorder="little")

        # Count trailing zeros
        trailing_zeros = (pow & -pow).bit_length() - 1

        if trailing_zeros >= proof_of_work_bits:
            sponge.absorb(nonce_bytes)
            sponge.squeeze_bytes(4)
            return nonce

        nonce += 1


def proof_of_work_verify(
    sponge, proof_of_work_bits: int, pow_nonce: int | None
) -> bool:
    assert proof_of_work_bits <= 32, "proof_of_work_bits must be <= 32"

    if proof_of_work_bits == 0:
        return True

    if pow_nonce is None:
        return False

    nonce_bytes = pow_nonce.to_bytes(length=8, byteorder="little")
    sponge.absorb(nonce_bytes)
    pow_bytes = sponge.squeeze_bytes(4)
    pow = int.from_bytes(pow_bytes, byteorder="little")

    # Count trailing zeros
    trailing_zeros = (pow & -pow).bit_length() - 1

    return trailing_zeros >= proof_of_work_bits
