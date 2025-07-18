import argparse
import sys
import os
from pathlib import Path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
from load_config_paths import PipelinePaths

from rdkit import Chem
import shutil

def merge_sdf_files(input_files, output_file):
    writer = Chem.SDWriter(str(output_file))
    for file in input_files:
        if not  os.path.exists(str(file)):
            print(f"Warning: File {file} does not exist. Skipping...")
            continue
        supp = Chem.SDMolSupplier(str(file))
        for mol in supp:
            if mol is not None:
                writer.write(mol)
    writer.close()




if __name__ == '__main__':
    paths = PipelinePaths()

    parser = argparse.ArgumentParser(description='Wrapper for CADD pipeline targeting Aurora protein kinases.')
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=False, default='B', help='Aurora kinase type (str, A, B)')
    parser.add_argument('--pdbid', type=str, required=False, default='4af3', help='Aurora kinase type (str, A, B)')
    parser.add_argument('--experiment', type=str, required=False, default='default', help='Aurora kinase type (str, A, B)')
    parser.add_argument('--output_file', type=str, required=False, default=None, help='Output file path')
    args = parser.parse_args()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    aurora = args.aurora
    pdbid = args.pdbid.lower()
    experiment = args.experiment

    # output_csv = paths.output_path(epoch, num_gen, known_binding_site, pdbid, args.output_file) 
    ligands_dir = Path(paths.equibind_ligands_path(experiment, epoch, num_gen, known_binding_site, pdbid))
    ligands_dir.mkdir(parents=True, exist_ok=True)

    output_sdf = ligands_dir / 'multiligand.sdf'
    print("Output sdf: ", output_sdf)

    # Get all SDF files in the ligands directory
    sdf_files = list(ligands_dir.glob("*.sdf"))
    
    if not sdf_files:
        print(f"Warning: No SDF files found in {ligands_dir}")
        # Create empty SDF file to satisfy Snakemake
        with open(output_sdf, 'w') as f:
            f.write("")
    else:
        print(f"Found {len(sdf_files)} SDF files to merge")
        merge_sdf_files(sdf_files, output_sdf)

    print(f"Merged SDF file created: {output_sdf}")

    # Copy multiligand.sdf to results directory using the new function
    results_sdf_path = paths.equibind_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'multiligand.sdf')

    try:
        shutil.copy2(output_sdf, results_sdf_path)
        print(f"Copied {output_sdf} to {results_sdf_path}")
    except Exception as e:
        print(f"Error copying file: {e}")
        # If copy fails, try creating the file directly in the results directory
        print("Attempting to create file directly in results directory...")
        try:
            if sdf_files:
                merge_sdf_files(sdf_files, results_sdf_path)
            else:
                with open(results_sdf_path, 'w') as f:
                    f.write("")
            print(f"Successfully created file directly: {results_sdf_path}")
        except Exception as e2:
            print(f"Failed to create file directly: {e2}")
            raise
