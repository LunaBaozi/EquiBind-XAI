#!/usr/bin/env python3
"""
EquiBind to Delta LinF9 Processor
Processes EquiBind output ligands through Delta LinF9 XGB scoring and ranks them.
"""

import os
import sys
import subprocess
import pandas as pd
import argparse
import time
from pathlib import Path
import tempfile
import re
from rdkit import Chem

def log(level, message):
    """Simple logging function"""
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    colors = {
        'INFO': '\033[0;32m',
        'WARN': '\033[1;33m', 
        'ERROR': '\033[0;31m',
        'DEBUG': '\033[0;34m'
    }
    reset = '\033[0m'
    color = colors.get(level, '')
    print(f"{color}[{level}]{reset} {timestamp} - {message}")

# Global variable for debug mode
DEBUG_MODE = False

def read_success_file(success_file):
    """Read the success.txt file and return list of successful ligand IDs"""
    ligand_ids = []
    try:
        with open(success_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.isdigit():
                    ligand_ids.append(int(line))
        return ligand_ids
    except Exception as e:
        log("ERROR", f"Failed to read success file: {e}")
        return []

def split_multiligand_sdf(sdf_file, output_dir, success_ids=None):
    """Split a multi-ligand SDF file into individual ligand files"""
    
    try:
        # Read the multi-ligand SDF file
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        ligand_files = []
        
        for i, mol in enumerate(suppl):
            if mol is None:
                log("WARN", f"Could not read molecule at index {i}")
                continue
                
            # If success_ids is provided, only process successful ligands
            if success_ids is not None and i not in success_ids:
                continue
            
            # Create individual SDF file for this ligand
            ligand_file = os.path.join(output_dir, f"ligand_{i}.sdf")
            writer = Chem.SDWriter(ligand_file)
            writer.write(mol)
            writer.close()
            
            ligand_files.append((i, ligand_file))
            if DEBUG_MODE:
                log("DEBUG", f"Extracted ligand {i} to {ligand_file}")
        
        log("INFO", f"Split multi-ligand SDF into {len(ligand_files)} individual files")
        return ligand_files
        
    except Exception as e:
        log("ERROR", f"Failed to split multi-ligand SDF: {e}")
        return []

def run_delta_linf9(protein_file, ligand_file, timeout=300):
    """Run Delta LinF9 on a single ligand"""
    
    # Path to Delta LinF9
    delta_linf9_dir = "/vol/data/drug-design-pipeline/external/deltalinf9"
    script_path = os.path.join(delta_linf9_dir, "script", "runXGB.py")
    
    # Convert to absolute paths
    protein_file = os.path.abspath(protein_file)
    ligand_file = os.path.abspath(ligand_file)
    
    try:
        # Change to Delta LinF9 directory
        original_dir = os.getcwd()
        os.chdir(delta_linf9_dir)
        
        # Run the script
        cmd = ["python", script_path, protein_file, ligand_file]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        # Change back to original directory
        os.chdir(original_dir)
        
        if result.returncode == 0:
            # Parse XGB score from output
            xgb_match = re.search(r'XGB \(in pK\) :\s*([0-9.-]+)', result.stdout)
            if xgb_match:
                xgb_score = float(xgb_match.group(1))
                return xgb_score, result.stdout
            else:
                log("WARN", "Could not parse XGB score from output")
                return None, result.stdout
        else:
            log("ERROR", f"Delta LinF9 failed: {result.stderr}")
            return None, result.stderr
            
    except subprocess.TimeoutExpired:
        log("ERROR", f"Delta LinF9 timed out after {timeout} seconds")
        return None, "TIMEOUT"
    except Exception as e:
        log("ERROR", f"Delta LinF9 execution failed: {e}")
        return None, str(e)

def process_equibind_output(equibind_dir, protein_file, output_dir, max_ligands=None):
    """Process all EquiBind output ligands through Delta LinF9"""
    
    ligands_dir = os.path.join(equibind_dir, "ligands")
    success_file = os.path.join(ligands_dir, "success.txt")
    multiligand_sdf = os.path.join(ligands_dir, "output.sdf")
    
    # Validate inputs
    if not os.path.exists(ligands_dir):
        log("ERROR", f"Ligands directory not found: {ligands_dir}")
        return False
        
    if not os.path.exists(multiligand_sdf):
        log("ERROR", f"Multi-ligand SDF file not found: {multiligand_sdf}")
        return False
        
    if not os.path.exists(protein_file):
        log("ERROR", f"Protein file not found: {protein_file}")
        return False
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create temporary directory for individual ligand files
    temp_dir = os.path.join(output_dir, "temp_ligands")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Read successful ligand IDs (if success file exists)
    success_ids = None
    if os.path.exists(success_file):
        success_ids = read_success_file(success_file)
        log("INFO", f"Found {len(success_ids)} successful ligands in success.txt")
    else:
        log("WARN", "No success.txt file found, processing all ligands in output.sdf")
    
    # Split the multi-ligand SDF file
    ligand_files = split_multiligand_sdf(multiligand_sdf, temp_dir, success_ids)
    
    if not ligand_files:
        log("ERROR", "No ligands found in the multi-ligand SDF file")
        return False
    
    # Limit number of ligands if specified
    if max_ligands:
        ligand_files = ligand_files[:max_ligands]
    
    log("INFO", f"Processing {len(ligand_files)} ligands through Delta LinF9")
    
    # Initialize results
    results = []
    
    # Process each ligand
    for i, (ligand_id, ligand_file) in enumerate(ligand_files):
        log("INFO", f"Processing ligand {ligand_id} ({i+1}/{len(ligand_files)})")
        
        start_time = time.time()
        xgb_score, output = run_delta_linf9(protein_file, ligand_file)
        processing_time = time.time() - start_time
        
        if xgb_score is not None:
            log("INFO", f"Ligand {ligand_id}: XGB score = {xgb_score:.3f}")
            results.append({
                'ligand_id': ligand_id,
                'ligand_file': ligand_file,
                'xgb_score': xgb_score,
                'processing_time': processing_time,
                'status': 'SUCCESS'
            })
        else:
            log("WARN", f"Failed to process ligand {ligand_id}")
            results.append({
                'ligand_id': ligand_id,
                'ligand_file': ligand_file,
                'xgb_score': None,
                'processing_time': processing_time,
                'status': 'FAILED'
            })
    
    # Create DataFrame and save results
    df = pd.DataFrame(results)
    
    # Save all results
    results_file = os.path.join(output_dir, "delta_linf9_results.csv")
    df.to_csv(results_file, index=False)
    log("INFO", f"Results saved to: {results_file}")
    
    # Create ranked results (only successful ones)
    successful_df = df[df['status'] == 'SUCCESS'].copy()
    ranked_file = os.path.join(output_dir, "ligands_ranked.csv")
    
    if len(successful_df) > 0:
        successful_df = successful_df.sort_values('xgb_score', ascending=False)
        successful_df.to_csv(ranked_file, index=False)
        log("INFO", f"Ranked results saved to: {ranked_file}")
        
        # Print top 10 results
        print(f"\nTOP 10 LIGANDS BY XGB SCORE:")
        print("-" * 60)
        print(f"{'Rank':>4} | {'Ligand ID':>9} | {'XGB Score':>10} | {'Processing Time':>15}")
        print("-" * 60)
        for idx, (_, row) in enumerate(successful_df.head(10).iterrows()):
            print(f"{idx+1:4} | {row['ligand_id']:9} | {row['xgb_score']:10.3f} | {row['processing_time']:13.2f}s")
    else:
        # Create empty ranked file if no successful ligands
        empty_df = pd.DataFrame(columns=['ligand_id', 'xgb_score', 'processing_time', 'status'])
        empty_df.to_csv(ranked_file, index=False)
        log("WARN", "No ligands were successfully processed")
        log("INFO", f"Empty ranked results saved to: {ranked_file}")
    
    # Clean up temporary files
    try:
        import shutil
        shutil.rmtree(temp_dir)
        log("INFO", "Cleaned up temporary ligand files")
    except Exception as e:
        log("WARN", f"Could not clean up temporary files: {e}")
    
    return True

def main():
    global DEBUG_MODE
    
    parser = argparse.ArgumentParser(description='Process EquiBind output through Delta LinF9')
    parser.add_argument('-i', '--input', required=True, 
                       help='EquiBind output directory (containing ligands/ subdirectory)')
    parser.add_argument('-p', '--protein', required=True,
                       help='Protein PDB file for scoring')
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory for results')
    parser.add_argument('-n', '--max-ligands', type=int,
                       help='Maximum number of ligands to process (for testing)')
    parser.add_argument('-d', '--debug', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Set debug mode
    DEBUG_MODE = args.debug
    
    log("INFO", "Starting EquiBind to Delta LinF9 processing")
    log("INFO", f"Input directory: {args.input}")
    log("INFO", f"Protein file: {args.protein}")
    log("INFO", f"Output directory: {args.output}")
    if args.max_ligands:
        log("INFO", f"Maximum ligands to process: {args.max_ligands}")
    
    success = process_equibind_output(args.input, args.protein, args.output, args.max_ligands)
    
    if success:
        log("INFO", "Processing completed successfully!")
        sys.exit(0)
    else:
        log("ERROR", "Processing failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
