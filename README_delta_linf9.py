#!/usr/bin/env python3
"""
EquiBind Delta LinF9 Pipeline Summary

Usage Examples:
1. Process all ligands from EquiBind output:
   python process_equibind_delta.py -i /path/to/equibind/output -p protein.pdb -o results

2. Process only first 5 ligands (for testing):
   python process_equibind_delta.py -i /path/to/equibind/output -p protein.pdb -o results -n 5

3. Process with debug output:
   python process_equibind_delta.py -i /path/to/equibind/output -p protein.pdb -o results -d

The script expects:
- EquiBind output directory containing ligands/output.sdf (multi-ligand file)
- Optionally ligands/success.txt (to filter successful ligands)
- A protein PDB file for Delta LinF9 scoring

Output files:
- delta_linf9_results.csv: All results
- ligands_ranked.csv: Results ranked by XGB score (best first)
"""

print(__doc__)
