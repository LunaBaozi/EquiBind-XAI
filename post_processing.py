import os, argparse
import pandas as pd
import shutil

from rdkit import Chem

from pathlib import Path
import sys

# Add project root to path to import from scripts
project_root = Path(__file__).parent.parent.parent
sys.path.append(str(project_root))

# try:
#     from scripts.aurk_int_preprocess import read_aurora_kinase_interactions
#     from scripts.gen_mols_preprocess import load_mols_from_sdf_folder
#     from scripts.load_config_paths import PipelinePaths
# except ImportError:
#     # Fallback imports if scripts not available
from load_config_paths import PipelinePaths

def load_confidence_scores(confidence_path):
    """
    Load confidence scores from EquiBind results.
    """
    if not os.path.exists(confidence_path):
        print(f"Warning: Confidence scores file not found at {confidence_path}")
        return pd.DataFrame()
    
    confidence_df = pd.read_csv(confidence_path)
    
    # The CSV now already has the filename column, no need to create it
    # Just ensure we have the expected columns
    if 'filename' not in confidence_df.columns:
        print("Warning: confidence_scores.csv missing 'filename' column")
        return pd.DataFrame()
    
    return confidence_df

def get_top_15_confidence_ligands(confidence_df):
    """
    Get the top 15 ligands by confidence score (higher is better).
    """
    if confidence_df.empty:
        print("Warning: No confidence scores available")
        return pd.DataFrame()
    
    # Sort by confidence score in descending order (higher confidence is better)
    top_15 = confidence_df.sort_values(by='confidence_score', ascending=False).head(15)
    
    return top_15

def merge_confidence_with_synthesizability(top_15_confidence, synth_path):
    """
    Merge top 15 confidence ligands with synthesizability scores.
    """
    if top_15_confidence.empty:
        print("Warning: No top confidence ligands available")
        return pd.DataFrame()
    
    if not os.path.exists(synth_path):
        print(f"Warning: Synthesizability scores file not found at {synth_path}")
        return top_15_confidence
    
    synth_df = pd.read_csv(synth_path)
    
    # Merge on filename
    merged_df = pd.merge(top_15_confidence, synth_df, on='filename', how='left')
    
    # Remove duplicates by keeping only the first occurrence of each filename
    # (This handles the case where hope-box has multiple Tanimoto scores per ligand)
    merged_df = merged_df.drop_duplicates(subset=['filename'], keep='first')
    
    # Sort by confidence score in descending order
    merged_df = merged_df.sort_values(by='confidence_score', ascending=False)
    
    # Ensure we only have 15 ligands
    merged_df = merged_df.head(15)
    
    return merged_df

def copy_top_15_sdf_files(top_15_df, source_sdf_path, dest_dir):
    """
    Copy the SDF files for the top 15 ligands to the destination directory.
    """
    os.makedirs(dest_dir, exist_ok=True)
    
    copied_files = 0
    copied_filenames = set()  # Track which filenames we've already copied
    
    # Read all molecules from output.sdf once
    source_file = Path(source_sdf_path).parent / 'output.sdf'
    
    if not os.path.exists(source_file):
        print(f"Error: Source SDF file not found: {source_file}")
        return 0
    
    try:
        # Load all molecules with their names
        supplier = Chem.SDMolSupplier(str(source_file))
        mol_dict = {}
        
        for mol in supplier:
            if mol is not None and mol.HasProp('_Name'):
                mol_name = mol.GetProp('_Name')
                mol_dict[mol_name] = mol
        
        print(f"Loaded {len(mol_dict)} molecules from {source_file}")
        
        # Copy molecules based on top 15 list
        for _, row in top_15_df.iterrows():
            filename = row['filename']
            
            # Skip if we've already copied this filename
            if filename in copied_filenames:
                continue
            
            # Use filename directly as it matches the _Name property in the SDF
            mol_name = filename
            
            if mol_name in mol_dict:
                mol = mol_dict[mol_name]
                dest_file = os.path.join(dest_dir, filename)
                
                # Write individual molecule to SDF file
                with Chem.SDWriter(dest_file) as writer:
                    writer.write(mol)
                
                copied_files += 1
                copied_filenames.add(filename)
                print(f"Copied molecule '{mol_name}' to {dest_file}")
            else:
                print(f"Warning: Molecule '{mol_name}' not found in output.sdf")
        
    except Exception as e:
        print(f"Error processing {source_file}: {e}")
        return 0
    
    print(f"Successfully copied {copied_files} SDF files to {dest_dir}")
    return copied_files

if __name__ == '__main__':
    paths = PipelinePaths()

    parser = argparse.ArgumentParser(description='Extract top 15 ligands by confidence score and prepare for Vina docking.')
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=False, default='B', help='Aurora kinase type (str, A, B)')
    parser.add_argument('--pdbid', type=str, required=False, default='4af3', help='PDB ID (str)')
    parser.add_argument('--experiment', type=str, required=False, default='default', help='Experiment name (str)')
    parser.add_argument('--output_file', type=str, required=False, default=None, help='Output file path')        
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    aurora = args.aurora
    pdbid = args.pdbid.lower()
    experiment = args.experiment

    # Step 1: Load confidence scores from EquiBind results
    confidence_path = paths.equibind_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'confidence_scores.csv')
    print(f"Loading confidence scores from: {confidence_path}")
    
    confidence_df = load_confidence_scores(confidence_path)
    
    if confidence_df.empty:
        print("No confidence scores found. Exiting.")
        exit(1)
    
    # Step 2: Get top 15 ligands by confidence score
    top_15_confidence = get_top_15_confidence_ligands(confidence_df)
    print(f"Found top 15 ligands by confidence score")
    
    # Step 3: Load synthesizability scores and merge
    synth_path = paths.hope_box_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'merged_scores.csv')
    print(f"Loading synthesizability scores from: {synth_path}")
    
    merged_results = merge_confidence_with_synthesizability(top_15_confidence, synth_path)
    
    # Step 4: Save combined results
    output_csv = paths.equibind_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'top_15_confidence_with_synth.csv')
    
    if not merged_results.empty:
        merged_results.to_csv(output_csv, index=False)
        print(f"Top 15 confidence ligands with synthesizability scores saved to: {output_csv}")
        
        # Display summary
        print(f"\nTop 15 ligands by confidence score:")
        print(merged_results[['filename', 'confidence_score', 'SA_score', 'tanimoto']].head(15))
    else:
        print("No merged results to save.")
        exit(1)
    
    # Step 4.5: Copy the CSV file to Vina-box directory
    vina_base_dir = os.path.join(paths.project_root, "external", "vina-box", pdbid, f"experiment_{experiment}_{epoch}_{num_gen}_{known_binding_site}_{pdbid}")
    vina_ligands_dir = os.path.join(vina_base_dir, "ligands")
    os.makedirs(vina_ligands_dir, exist_ok=True)
    
    # Copy the CSV file to the base experiment directory (not ligands subdirectory)
    vina_csv_dest = os.path.join(vina_base_dir, 'top_15_confidence_with_synth.csv')
    shutil.copy2(output_csv, vina_csv_dest)
    print(f"Copied CSV file to Vina directory: {vina_csv_dest}")

    # Step 4.6: Copy the PDB file to the Vina-box base directory for the protein
    vina_protein_dir = os.path.join(paths.project_root, "external", "vina-box", pdbid)
    os.makedirs(vina_protein_dir, exist_ok=True)

    # Find the PDB file dynamically based on current experiment
    pdb_source = paths.equibind_data_path(experiment, epoch, num_gen, known_binding_site, pdbid)
    pdb_source = os.path.join(pdb_source, "ligands", f"{pdbid}_A_rec_reduce_noflip.pdb")
    pdb_dest = os.path.join(vina_protein_dir, f"{pdbid}_A_rec_reduce_noflip.pdb")

    try:
        if os.path.exists(pdb_source):
            shutil.copy2(pdb_source, pdb_dest)
            print(f"Copied PDB file to Vina directory: {pdb_dest}")
        else:
            print(f"Warning: PDB source file not found: {pdb_source}")
            # Fallback: try to find it in any experiment directory
            search_pattern = os.path.join(paths.project_root, f"external/equibind/data/{pdbid}/*/ligands/{pdbid}_A_rec_reduce_noflip.pdb")
            import glob
            matching_files = glob.glob(search_pattern)
            if matching_files:
                pdb_source = matching_files[0]  # Use the first match
                shutil.copy2(pdb_source, pdb_dest)
                print(f"Found and copied PDB file from: {pdb_source} to {pdb_dest}")
            else:
                print(f"Error: Could not find PDB file for {pdbid}")
    except Exception as e:
        print(f"Warning: Could not copy PDB file: {e}")


    # Step 5: Copy SDF files to Vina-box directory
    # Create vina-box destination directory (already created above)
    # vina_dest_dir = f"/vol/data/drug-design-pipeline/external/vina-box/docking/{pdbid}/experiment_{experiment}_{epoch}_{num_gen}_{known_binding_site}_{pdbid}"
    # Source SDF path (from EquiBind results)
    source_sdf_path = paths.equibind_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'output.sdf')
    
    print(f"\nCopying SDF files to Vina ligands directory: {vina_ligands_dir}")
    copied_count = copy_top_15_sdf_files(merged_results, source_sdf_path, vina_ligands_dir)

    if copied_count > 0:
        print(f"Successfully prepared {copied_count} ligands for Vina docking")
    else:
        print("Warning: No SDF files were copied")
    
    print(f"\nProcessing complete!")
    print(f"Results CSV: {output_csv}")
    print(f"Vina CSV: {vina_csv_dest}")
    print(f"SDF files location: {vina_ligands_dir}")
    print(f"Total ligands processed: {len(merged_results)}")


