from dataclasses import dataclass
from enum import Enum


class SoundnessType(Enum):
    PROVABLE = "provable"
    CONJECTURE = "conjecture"


@dataclass
class Parameters:
    security_level: int
    protocol_security_level: int
    starting_degree: int
    stopping_degree: int
    folding_factor: int
    starting_rate: int
    soundness_type: SoundnessType


@dataclass
class FullParameters(Parameters):
    num_rounds: int
    rates: list
    repetitions: list[int]
    pow_bits: list[int]
    ood_samples: int
    degrees: list


DEFAULT_PARAMETERS = Parameters(
    starting_degree=2**10,
    stopping_degree=2**6,
    protocol_security_level=106,
    security_level=128,
    starting_rate=2,
    folding_factor=16,
    soundness_type=SoundnessType.CONJECTURE,
)
