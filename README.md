# Hash-based Proof-carrying Data Schemes (WIP)

This repo is work-in-progress implementation of several hash-based proof-carrying data schemes. 

Proof-carrying data (PCD) is a cryptographic primitive for computational integrity in a distributed setting. State-of-the-art constructions of PCD are based on accumulation schemes (and, closely related, folding schemes).

## Implemented details

The following constructions will be implemented:

- **STIR**, ([paper](https://eprint.iacr.org/2024/390)). While STIR is not PCD, its ideas have been used in to build PCD. This implementation will be a re-write of [Rust implmentation](https://github.com/WizardOfMenlo/stir) by the autors.
- **ARC**, ([paper](https://eprint.iacr.org/2024/1731)). ARC is hash-based accumulation scheme that supports an unbounded number of accumulation steps.
- **WARP**, ([paper](https://eprint.iacr.org/2025/753)). WARP, the first accumulation scheme with linear prover time and logarithmic verifier time.

WARNING: This is intended to be an academic prototype and by no means is suitable for production use.

