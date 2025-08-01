import math

n = 54


def grain_sr_generator():
    bit_sequence = [1] * 80
    for _ in range(0, 160):
        new_bit = (
            bit_sequence[62]
            ^ bit_sequence[51]
            ^ bit_sequence[38]
            ^ bit_sequence[23]
            ^ bit_sequence[13]
            ^ bit_sequence[0]
        )
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)

    while True:
        new_bit = (
            bit_sequence[62]
            ^ bit_sequence[51]
            ^ bit_sequence[38]
            ^ bit_sequence[23]
            ^ bit_sequence[13]
            ^ bit_sequence[0]
        )
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        while new_bit == 0:
            new_bit = (
                bit_sequence[62]
                ^ bit_sequence[51]
                ^ bit_sequence[38]
                ^ bit_sequence[23]
                ^ bit_sequence[13]
                ^ bit_sequence[0]
            )
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
            new_bit = (
                bit_sequence[62]
                ^ bit_sequence[51]
                ^ bit_sequence[38]
                ^ bit_sequence[23]
                ^ bit_sequence[13]
                ^ bit_sequence[0]
            )
            bit_sequence.pop(0)
            bit_sequence.append(new_bit)
        new_bit = (
            bit_sequence[62]
            ^ bit_sequence[51]
            ^ bit_sequence[38]
            ^ bit_sequence[23]
            ^ bit_sequence[13]
            ^ bit_sequence[0]
        )
        bit_sequence.pop(0)
        bit_sequence.append(new_bit)
        yield new_bit


grain_gen = grain_sr_generator()


def grain_random_bits(num_bits):
    random_bits = [next(grain_gen) for i in range(0, num_bits)]
    random_int = int("".join(str(i) for i in random_bits), 2)
    return random_int


def gen_constants():
    # Round constants for GF(2^n)
    # num_rounds = math.ceil(n // math.log(3, 2))
    num_rounds = 54
    round_constants = [0x0]
    for i in range(1, num_rounds):
        random_int = grain_random_bits(n)
        round_constants.append(random_int)
    print("len(round_constants)", len(round_constants))
    print("Round constants for GF(2^n):")
    hex_length = int(math.ceil(float(n) / 4)) + 2  # +2 for "0x"
    print(["{0:#0{1}x}".format(entry, hex_length) for entry in round_constants])
    for c in round_constants:
        assert c < 2**n - 1
