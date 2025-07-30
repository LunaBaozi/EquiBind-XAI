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
    
    for _, row in top_15_df.iterrows():
        filename = row['filename']
        
        # Skip if we've already copied this filename
        if filename in copied_filenames:
            continue
        
        # Extract ligand index from filename (e.g., "0" from "0.sdf")
        ligand_index = int(filename.replace('.sdf', ''))
        
        # Try to find the source SDF file
        # First try the output.sdf file and extract individual molecules
        source_file = Path(source_sdf_path).parent / 'output.sdf'
        
        if os.path.exists(source_file):
            try:
                # Read all molecules from output.sdf
                supplier = Chem.SDMolSupplier(str(source_file))
                mols = [mol for mol in supplier if mol is not None]
                
                if ligand_index < len(mols):
                    mol = mols[ligand_index]
                    dest_file = os.path.join(dest_dir, filename)
                    
                    # Write individual molecule to SDF file
                    with Chem.SDWriter(dest_file) as writer:
                        writer.write(mol)
                    
                    copied_files += 1
                    copied_filenames.add(filename)
                    print(f"Copied ligand {ligand_index} ({filename}) to {dest_file}")
                else:
                    print(f"Warning: Ligand index {ligand_index} not found in output.sdf")
            
            except Exception as e:
                print(f"Error copying ligand {filename}: {e}")
        else:
            print(f"Warning: Source SDF file not found: {source_file}")
    
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
    
    # Step 5: Copy SDF files to Vina-box directory
    # Create vina-box destination directory
    vina_dest_dir = f"/vol/data/drug-design-pipeline/external/vina-box/docking/{pdbid}/experiment_{experiment}_{epoch}_{num_gen}_{known_binding_site}_{pdbid}"
    
    # Source SDF path (from EquiBind results)
    source_sdf_path = paths.equibind_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'output.sdf')
    
    print(f"\nCopying SDF files to Vina directory: {vina_dest_dir}")
    copied_count = copy_top_15_sdf_files(merged_results, source_sdf_path, vina_dest_dir)
    
    if copied_count > 0:
        print(f"Successfully prepared {copied_count} ligands for Vina docking")
    else:
        print("Warning: No SDF files were copied")
    
    print(f"\nProcessing complete!")
    print(f"Results CSV: {output_csv}")
    print(f"SDF files location: {vina_dest_dir}")
    print(f"Total ligands processed: {len(merged_results)}")


