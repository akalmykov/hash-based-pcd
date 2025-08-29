import sys
import os

# Ensure the current directory is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from combine import test_combine
from rs_accumulation import test_rs_accumulation
from r1cs_to_racc import test_r1cs_to_racc

if __name__ == "__main__":
    test_combine()
    test_r1cs_to_racc()
    test_rs_accumulation()
