import sys
import os

# Ensure the current directory is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from combine import test_combine
from rs_accumulation import test_rs_accumulation

if __name__ == "__main__":
    # test_combine()
    test_rs_accumulation()
