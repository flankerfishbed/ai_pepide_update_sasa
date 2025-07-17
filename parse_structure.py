from typing import Dict, List, Any
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1

def parse_pdb_structure(pdb_content: str, chain_id: str = 'A') -> Dict[str, Any]:
    """
    Parse PDB content and extract:
      - Primary sequence
      - Residue list with CA atom coordinates
      - Residue type, ID, and 3D position
    Only standard amino acids are included.
    Returns a dictionary with extracted data.
    """
    import io
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', io.StringIO(pdb_content))
    model = structure[0]
    if chain_id not in model:
        raise ValueError(f"Chain {chain_id} not found in structure.")
    chain = model[chain_id]

    residues = []
    sequence = ''
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resname = res.get_resname().capitalize()
        res_id = None
        res_id_tuple = res.get_id()
        if res_id_tuple is not None and isinstance(res_id_tuple, (tuple, list)) and len(res_id_tuple) > 1:
            possible_id = res_id_tuple[1]
            if isinstance(possible_id, int):
                res_id = possible_id
        ca_atom = res['CA'] if 'CA' in res else None
        coords = ca_atom.get_coord().tolist() if ca_atom is not None else None
        residues.append({
            'name': resname,
            'id': res_id,
            'ca_coord': coords
        })
        # Use mapping for three-letter to one-letter code
        one_letter = protein_letters_3to1.get(resname, None)
        if one_letter:
            sequence += one_letter

    return {
        'sequence': sequence,
        'residues': residues
    } 