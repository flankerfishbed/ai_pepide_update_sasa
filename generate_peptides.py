import random
from typing import List, Optional

# Standard amino acids
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

def generate_random_peptides(
    sequence: str,
    num_peptides: int = 5,
    peptide_length: int = 8,
    charge: Optional[str] = None,
    hydrophobicity: Optional[str] = None
) -> List[str]:
    """
    Generate a list of random peptide sequences from the given protein sequence.
    Optionally filter by charge or hydrophobicity (stub for now).
    """
    if not sequence or num_peptides <= 0 or peptide_length <= 0:
        return []
    peptides = set()
    max_start = max(0, len(sequence) - peptide_length)
    while len(peptides) < num_peptides and max_start > 0:
        start = random.randint(0, max_start)
        peptide = sequence[start:start+peptide_length]
        if len(peptide) == peptide_length:
            peptides.add(peptide)
    return list(peptides) 