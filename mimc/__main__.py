import argparse

from .mimc_prime_field import mimc_main

from .mimc_constants_generation import gen_constants


def main():
    mimc_main()
    # gen_constants()


if __name__ == "__main__":
    main()
