import freesasa
import pandas as pd
from typing import Dict, Any, List, Tuple

def analyze_surface_residues(pdb_content: str, selected_chain: str = 'A') -> pd.DataFrame:
    """
    Placeholder for surface residue analysis (SASA-based analysis removed).
    Returns all residues in the selected chain as 'surface' residues with dummy SASA values.
    
    Args:
        pdb_content (str): PDB file content as string
        selected_chain (str): Chain ID to analyze (default: 'A')
    
    Returns:
        pd.DataFrame: DataFrame with columns ['Residue Name', 'Residue Number', 'Chain', 'SASA', 'Property']
    """
    # Simple PDB parsing to extract residues for the selected chain
    residues = set()
    for line in pdb_content.splitlines():
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chain = line[21].strip()
            if chain != selected_chain:
                continue
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            residues.add((resname, resnum, chain))
    
    # Dummy classification and SASA
    def classify_residue(residue_name: str) -> str:
        HYDROPHOBIC = {'A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W'}
        CHARGED = {'D', 'E', 'K', 'R', 'H'}
        if residue_name[0] in HYDROPHOBIC:
            return "Hydrophobic"
        elif residue_name[0] in CHARGED:
            return "Charged"
        else:
            return "Polar/Other"
    
    residues_data = []
    for resname, resnum, chain in residues:
        property_type = classify_residue(resname)
        residues_data.append({
            'Residue Name': resname,
            'Residue Number': resnum,
            'Chain': chain,
            'SASA': 50.0,  # Dummy value
            'Property': property_type
        })
    df = pd.DataFrame(residues_data)
    if not df.empty:
        df = df.sort_values('Residue Number')
    return df

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
    Args:
        pdb_filepath (str): Path to the PDB file.
        selected_chain (str): Chain ID to analyze (default: 'A').
        sasa_threshold (float): Minimum SASA value to consider a residue exposed.
    Returns:
        List[Dict[str, Any]]: List of dicts with keys: 'chain', 'residue_number', 'residue_name', 'sasa'
    """
    structure = freesasa.Structure(pdb_filepath)
    result = freesasa.calc(structure)
    exposed_residues = []
    for key, sasa in result.residueAreas().items():
        # key is (chain, resnum, resname)
        if isinstance(key, tuple) and len(key) == 3:
            chain, resnum, resname = key
            try:
                resnum = int(resnum)
            except ValueError:
                continue  # skip if not a valid integer
        else:
            continue
        if chain != selected_chain:
            continue
        if sasa > sasa_threshold:
            exposed_residues.append({
                'chain': chain,
                'residue_number': resnum,  # int
                'residue_name': resname,
                'sasa': float(sasa)        # float
            })
    return exposed_residues 