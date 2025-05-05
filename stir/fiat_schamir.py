import blake3
import galois
from typing import List


class Blake3Sponge:
    def __init__(self):
        self.hasher = blake3.blake3()

    def absorb(self, input_data: bytes):
        self.hasher.update(input_data)

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

    def squeeze_field(self, GF, size):
        bits = self.squeeze_bits(size)
        value = 0
        for i, bit in enumerate(bits):
            if bit:
                value += 1 << i
        return GF(value)


# Example Usage:
if __name__ == "__main__":
    sponge = Blake3Sponge()

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

    with open("../../stir/fs1.json", "r") as f:
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
    sponge = Blake3Sponge()
    sponge.absorb(bytes(merkle_root))
    fs1 = sponge.squeeze_bytes(1024)
    fs1_list = list(fs1)
    assert fs1_list == fs

    with open("../../stir/fs_bits.json", "r") as f:
        t_fs_bits = json.load(f)
    sponge = Blake3Sponge()
    sponge.absorb(bytes(merkle_root))
    fs_bits = sponge.squeeze_bits(1024)
    assert t_fs_bits == fs_bits

    sponge = Blake3Sponge()
    sponge.absorb(bytes(merkle_root))
    from constants import FIELD192

    GF = galois.GF(FIELD192)

    folding_randomness = sponge.squeeze_field(GF, 191)
    with open("../../stir/folding_randomness.json", "r") as f:
        t_folding_randomness = int(json.load(f))

    assert folding_randomness == t_folding_randomness
