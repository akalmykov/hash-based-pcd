from stir.stir_main import run_stir
import argparse

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run STIR protocol with options")
    parser.add_argument(
        "--prove",
        dest="prove",
        action="store_true",
        default=True,
        help="Run the proving part of the protocol (default: True)",
    )
    parser.add_argument(
        "--no-prove",
        dest="prove",
        action="store_false",
        help="Skip the proving part of the protocol",
    )
    parser.add_argument(
        "--verify",
        dest="verify",
        action="store_true",
        default=True,
        help="Run the verification part of the protocol (default: True)",
    )
    parser.add_argument(
        "--no-verify",
        dest="verify",
        action="store_false",
        help="Skip the verification part of the protocol",
    )

    args = parser.parse_args()

    # Pass the parsed arguments to run_stir
    run_stir(prove=args.prove, verify=args.verify)
