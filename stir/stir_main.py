import galois
import time
from constants import FIELD192


def main():
    start_time = time.time()
    field = galois.GF(FIELD192)
    end_time = time.time()
    print(f"Generate FIELD192: {end_time - start_time} seconds")

    starting_degree = 2**20
    stopping_degree = 2**6
    protocol_security_level = 106
    security_level = 128
    starting_rate = 2
    folding_factor = 16

    start_time = time.time()
    poly = field.Random(starting_degree + 1, seed=1)
    end_time = time.time()
    print(f"Generate random polynomial: {end_time - start_time} seconds")


if __name__ == "__main__":
    main()
