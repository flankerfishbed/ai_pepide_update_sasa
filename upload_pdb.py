from typing import Optional
import os

def validate_pdb_file(filename: str) -> bool:
    """
    Check if the uploaded file is a valid .pdb file by extension.
    """
    return filename.lower().endswith('.pdb')

def read_pdb_file(file_path: str) -> Optional[str]:
    """
    Read the contents of a PDB file and return as a string.
    Returns None if file does not exist or is not readable.
    """
    if not os.path.isfile(file_path):
        return None
    try:
        with open(file_path, 'r') as f:
            return f.read()
    except Exception as e:
        print(f"Error reading file: {e}")
        return None 