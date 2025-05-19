import blake3
import galois
from typing import List
from .merkle import serialize_field_vector


class Blake3Sponge:

    def __init__(self, GF, field_size):
        self.hasher = blake3.blake3()
        self.GF = GF
        self.field_size = field_size

    def clone(self):
        sponge_clone = Blake3Sponge(self.GF, self.field_size)
        sponge_clone.hasher = self.hasher.copy()
        return sponge_clone

    def absorb(self, input_data: bytes):
        self.hasher.update(input_data)

    def absorb_field_vec(self, vector: list):
        bytes_to_absorb = serialize_field_vector(vector)[8:]
        # print("bytes_to_absorb", list(bytes_to_absorb[8:]))
        self.absorb(bytes_to_absorb)

    def squeeze_bytes(self, num_bytes: int) -> bytes:
        digest = self.hasher.digest(length=num_bytes)
        self.hasher.update(digest)
        return digest

    def squeeze_bits(self, num_bits: int) -> List[bool]:
        num_bytes = (num_bits + 7) // 8
        byte_output = self.squeeze_bytes(num_bytes)  # This already updates the state

        bits = []
        for byte in byte_output:
            for i in range(8):
                bits.append((byte >> (7 - i)) & 1 == 1)

        return bits[:num_bits]

    def reset(self):
        self.hasher = blake3.blake3()

    def squeeze_f(self):
        return self.squeeze_field_element(self.GF, self.field_size)

    def squeeze_int(self, mod):
        assert (mod & (mod - 1) == 0) and mod > 0, "Range must be a power of 2"

        bytes_array = self.squeeze_bytes(8)
        candidate = int.from_bytes(bytes_array, byteorder="little")
        return candidate % mod

    def squeeze_field_element(self, GF, size):
        bits = self.squeeze_bits(size)
        value = 0
        for i, bit in enumerate(bits):
            if bit:
                value += 1 << i
        return GF(value)

    def squeeze_several_f(self, n):
        return self.squeeze_field_elements(self.GF, self.field_size, n)

    def squeeze_field_elements(self, GF, size, n):
        elements = []
        bits = self.squeeze_bits(size * n)
        for j in range(n):
            value = 0
            for i, bit in enumerate(bits[size * j : size * (j + 1)]):
                if bit:
                    value += 1 << i
            elements.append(GF(value))
        return elements


# Example Usage:
if __name__ == "__main__":
    from constants import FIELD192

    GF = galois.GF(FIELD192)

    sponge = Blake3Sponge(GF, 191)

    # Absorb some data
    sponge.absorb(b"some initial data")
    sponge.absorb(b"more data")

    # Squeeze bytes
    squeezed_bytes = sponge.squeeze_bytes(3200)
    print(f"Squeezed bytes: {squeezed_bytes.hex()}")

    # Squeeze more bytes (state has been updated)
    more_squeezed_bytes = sponge.squeeze_bytes(16)
    print(f"More squeezed bytes: {more_squeezed_bytes.hex()}")

    # Reset the sponge
    sponge.reset()
    sponge.absorb(b"data for bit squeezing")

    # Squeeze bits
    squeezed_bits = sponge.squeeze_bits(100)
    print(
        f"Squeezed bits ({len(squeezed_bits)}): {''.join(['1' if b else '0' for b in squeezed_bits])}"
    )

    # Squeeze more bits (state has been updated)
    more_squeezed_bits = sponge.squeeze_bits(50)
    print(
        f"More squeezed bits ({len(more_squeezed_bits)}): {''.join(['1' if b else '0' for b in more_squeezed_bits])}"
    )
    import json

    with open("./test/original-stir-traces/fs1.json", "r") as f:
        fs = json.load(f)

    merkle_root = [
        174,
        31,
        159,
        28,
        72,
        126,
        159,
        221,
        10,
        84,
        150,
        145,
        57,
        183,
        55,
        78,
        213,
        251,
        157,
        12,
        184,
        208,
        155,
        246,
        208,
        235,
        64,
        206,
        58,
        128,
        85,
        247,
    ]

    sponge = Blake3Sponge(GF, 191)
    sponge.absorb(bytes(merkle_root))
    fs1 = sponge.squeeze_bytes(1024)
    fs1_list = list(fs1)
    assert fs1_list == fs

    with open("./test/original-stir-traces/fs_bits.json", "r") as f:
        t_fs_bits = json.load(f)
    sponge = Blake3Sponge(GF, 191)
    sponge.absorb(bytes(merkle_root))
    fs_bits = sponge.squeeze_bits(1024)
    assert t_fs_bits == fs_bits

    sponge = Blake3Sponge(GF, 191)
    sponge.absorb(bytes(merkle_root))

    folding_randomness = sponge.squeeze_field_element(GF, 191)
    with open("./test/original-stir-traces/folding_randomness.json", "r") as f:
        t_folding_randomness = int(json.load(f))

    assert folding_randomness == t_folding_randomness
