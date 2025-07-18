import freesasa
import pandas as pd
from typing import Dict, Any, List, Tuple

def analyze_surface_residues(pdb_content: str, selected_chain: str = 'A') -> pd.DataFrame:
    """
    Analyze surface residues using FreeSASA to calculate solvent-accessible surface area (SASA).
    Only standard amino acids are included.
    Returns a DataFrame with columns ['Residue Name', 'Residue Number', 'Chain', 'SASA', 'Property']
    """
    # Write PDB content to a temporary file
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as temp_file:
        temp_file.write(pdb_content)
        temp_file_path = temp_file.name

    # Standard amino acids (3-letter codes)
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    HYDROPHOBIC = {'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W'}
    CHARGED = {'D', 'E', 'K', 'R', 'H'}

    def classify_residue(residue_name: str) -> str:
        if residue_name[0] in HYDROPHOBIC:
            return "Hydrophobic"
        elif residue_name[0] in CHARGED:
            return "Charged"
        else:
            return "Polar/Other"

    try:
        structure = freesasa.Structure(temp_file_path)
        result = freesasa.calc(structure)
        residues_data = []
        print("--- DEBUG: All residues and SASA values (before filtering) ---")
        for key, sasa in result.residueAreas().items():
            print(f"key={key}, sasa={sasa}")
        print("--- END DEBUG ---")
        for key, sasa in result.residueAreas().items():
            if isinstance(key, tuple) and len(key) == 3:
                chain, resnum, resname = key
                try:
                    resnum = int(resnum)
                except ValueError:
                    continue
            else:
                continue
            if chain != selected_chain:
                continue
            if resname not in standard_residues:
                continue
            if sasa > 0.0:
                property_type = classify_residue(resname)
                residues_data.append({
                    'Residue Name': resname,
                    'Residue Number': resnum,
                    'Chain': chain,
                    'SASA': float(sasa),
                    'Property': property_type
                })
        df = pd.DataFrame(residues_data)
        if not df.empty:
            df = df.sort_values('Residue Number')
        return df
    finally:
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)

def get_surface_summary(surface_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Generate a summary of surface analysis results.
    
    Args:
        surface_df (pd.DataFrame): Surface analysis DataFrame
    
    Returns:
        Dict[str, Any]: Summary statistics
    """
    if surface_df.empty:
        return {
            'total_surface_residues': 0,
            'hydrophobic_count': 0,
            'charged_count': 0,
            'polar_count': 0,
            'avg_sasa': 0.0,
            'max_sasa': 0.0
        }
    summary = {
        'total_surface_residues': len(surface_df),
        'hydrophobic_count': len(surface_df[surface_df['Property'] == 'Hydrophobic']),
        'charged_count': len(surface_df[surface_df['Property'] == 'Charged']),
        'polar_count': len(surface_df[surface_df['Property'] == 'Polar/Other']),
        'avg_sasa': round(surface_df['SASA'].mean(), 2),
        'max_sasa': round(surface_df['SASA'].max(), 2)
    }
    return summary 

def get_exposed_residues_from_pdb(pdb_filepath: str, selected_chain: str = 'A', sasa_threshold: float = 25.0) -> List[Dict[str, Any]]:
    """
    Compute residue-level SASA using FreeSASA and return a list of exposed residues above a threshold.
    Only standard amino acids are included.
    """
    standard_residues = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    structure = freesasa.Structure(pdb_filepath)
    result = freesasa.calc(structure)
    exposed_residues = []
    for key, sasa in result.residueAreas().items():
        if isinstance(key, tuple) and len(key) == 3:
            chain, resnum, resname = key
            try:
                resnum = int(resnum)
            except ValueError:
                continue
        else:
            continue
        if chain != selected_chain:
            continue
        if resname not in standard_residues:
            continue
        if sasa > sasa_threshold:
            exposed_residues.append({
                'chain': chain,
                'residue_number': resnum,
                'residue_name': resname,
                'sasa': float(sasa)
            })
    return exposed_residues 
