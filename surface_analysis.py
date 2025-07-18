diff --git a/surface_analysis.py b/surface_analysis.py
index ef6ad0b2309e88e1d62e66970081fb5a85baf034..e59591474e48bbaae78b60c8a9e1af41ef59192e 100644
--- a/surface_analysis.py
+++ b/surface_analysis.py
@@ -13,69 +13,68 @@ def analyze_surface_residues(pdb_content: str, selected_chain: str = 'A') -> Tup
     with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as temp_file:
         temp_file.write(pdb_content)
         temp_file_path = temp_file.name
 
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
 
     debug_info = []
     try:
         structure = freesasa.Structure(temp_file_path)
         result = freesasa.calc(structure)
         residues_data = []
         residue_areas = result.residueAreas()
-        for chain_id, residues in residue_areas.items():
-            debug_info.append(f"chain={chain_id}, residues={list(residues.keys())}")
-            if chain_id != selected_chain:
-                continue
-            for resnum, area in residues.items():
+        chain_residues = residue_areas.get(selected_chain, {})
+        if not chain_residues:
+            debug_info.append(f"chain {selected_chain} not found")
+        for resnum, area in chain_residues.items():
                 try:
                     resnum_int = int(resnum)
                 except ValueError:
                     continue
                 resname = area.residueType
                 sasa = area.total
                 if resname not in standard_residues:
                     continue
                 if sasa > 0.0:
                     property_type = classify_residue(resname)
                     residues_data.append({
                         'Residue Name': resname,
                         'Residue Number': resnum_int,
-                        'Chain': chain_id,
+                        'Chain': selected_chain,
                         'SASA': float(sasa),
                         'Property': property_type
                     })
         df = pd.DataFrame(residues_data)
         if not df.empty:
             df = df.sort_values('Residue Number')
         return df, debug_info
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
diff --git a/surface_analysis.py b/surface_analysis.py
index ef6ad0b2309e88e1d62e66970081fb5a85baf034..e59591474e48bbaae78b60c8a9e1af41ef59192e 100644
--- a/surface_analysis.py
+++ b/surface_analysis.py
@@ -84,45 +83,43 @@ def get_surface_summary(surface_df: pd.DataFrame) -> Dict[str, Any]:
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
-    for chain_id, residues in result.residueAreas().items():
-        if chain_id != selected_chain:
-            continue
-        for resnum, area in residues.items():
+    chain_residues = result.residueAreas().get(selected_chain, {})
+    for resnum, area in chain_residues.items():
             try:
                 resnum_int = int(resnum)
             except ValueError:
                 continue
             resname = area.residueType
             sasa = area.total
             if resname not in standard_residues:
                 continue
             if sasa > sasa_threshold:
                 exposed_residues.append({
-                    'chain': chain_id,
+                    'chain': selected_chain,
                     'residue_number': resnum_int,
                     'residue_name': resname,
                     'sasa': float(sasa)
                 })
     return exposed_residues 
