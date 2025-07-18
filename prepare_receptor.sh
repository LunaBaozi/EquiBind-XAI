#!/bin/bash

# Usage: ./prepare_receptor.sh --epoch <epoch> --num_gen <num_gen> --known_binding_site <known_binding_site> --pdbid <pdbid> --aurora <aurora> --experiment <experiment>
set -e

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --epoch)
            epoch="$2"
            shift 2
            ;;
        --num_gen)
            num_gen="$2"
            shift 2
            ;;
        --known_binding_site)
            known_binding_site="$2"
            shift 2
            ;;
        --pdbid)
            pdbid="$2"
            shift 2
            ;;
        --aurora)
            aurora="$2"
            shift 2
            ;;
        --experiment)
            experiment="$2"
            shift 2
            ;;
        *)
            echo "Unknown option $1"
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$epoch" ] || [ -z "$num_gen" ] || [ -z "$known_binding_site" ] || [ -z "$pdbid" ] || [ -z "$aurora" ] || [ -z "$experiment" ]; then
    echo "Usage: $0 --epoch <epoch> --num_gen <num_gen> --known_binding_site <known_binding_site> --pdbid <pdbid> --aurora <aurora> --experiment <experiment>"
    exit 1
fi

src="data/${pdbid}/${pdbid}_A_rec.pdb"
ligands_dir="../../results/experiment_${experiment}_${epoch}_${num_gen}_${known_binding_site}_${pdbid}/ligands"
dst="${ligands_dir}/${pdbid}_A_rec.pdb"
output_file="${ligands_dir}/${pdbid}_A_rec_reduce_noflip.pdb"

echo "Looking for source file: $src"
if [ ! -f "$src" ]; then
    echo "Source file not found: $src"
    exit 2
fi

# Create the ligands directory if it doesn't exist
echo "Creating ligands directory: $ligands_dir"
mkdir -p "$ligands_dir"

# Copy (don't move) the file to preserve the original
echo "Copying $src to $dst"
cp "$src" "$dst"

# Check if reduce command is available
if ! command -v reduce &> /dev/null; then
    echo "WARNING: reduce command not found. Copying original file instead."
    cp "$dst" "$output_file"
else
    echo "Running reduce on $dst"
    echo "Output will be saved to: $output_file"
    
    # Check reduce version and help
    echo "Reduce version:"
    reduce --version 2>&1 || reduce -h 2>&1 | head -5
    
    # Debug: show current directory and file existence
    echo "Current directory: $(pwd)"
    echo "Input file exists: $(ls -la "$dst" 2>/dev/null || echo "NOT FOUND")"
    
    # Temporarily disable strict mode for the reduce command
    set +e
    reduce -NOFLIP "$dst" > "$output_file"
    reduce_exit_code=$?
    set -e
    
    if [ $reduce_exit_code -eq 0 ]; then
        echo "reduce command completed successfully"
    else
        echo "reduce -NOFLIP failed with exit code: $reduce_exit_code"
        echo "Trying without flags..."
        
        reduce "$dst" > "$output_file"
        reduce_exit_code=$?
        
        if [ $reduce_exit_code -eq 0 ]; then
            echo "reduce command completed successfully (without flags)"
        else
            echo "All reduce attempts failed, copying original file instead"
            cp "$dst" "$output_file"
        fi
    fi
fi

# Verify the output file was created
if [ -f "$output_file" ]; then
    echo "Output file created successfully: $output_file"
    ls -la "$output_file"
else
    echo "ERROR: Output file was not created: $output_file"
    exit 3
fi