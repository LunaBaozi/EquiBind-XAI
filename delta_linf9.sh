#!/bin/bash

# EquiBind to Delta LinF9 Processing Script
# This script takes EquiBind output ligands and processes them through the Delta LinF9 algorithm
# Then ranks ligands by their XGB binding affinity scores in descending order

set -e  # Exit on any error

# Default parameters
EQUIBIND_OUTPUT_DIR=""
PROTEIN_FILE=""
OUTPUT_DIR=""
CONDA_ENV="delta_linf9"
VERBOSE=false

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print usage
usage() {
    echo "Usage: $0 -i <equibind_output_dir> -p <protein_file> -o <output_dir> [options]"
    echo ""
    echo "Required arguments:"
    echo "  -i, --input <dir>      EquiBind output directory containing ligands/ subdirectory"
    echo "  -p, --protein <file>   Protein PDB file to use for scoring"
    echo "  -o, --output <dir>     Output directory for results"
    echo ""
    echo "Optional arguments:"
    echo "  -e, --env <name>       Conda environment name (default: delta_linf9)"
    echo "  -v, --verbose          Enable verbose output"
    echo "  -h, --help             Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -i /path/to/equibind/data/4af3/experiment_quasi_99_5_True_4af3 \\"
    echo "     -p /path/to/protein.pdb \\"
    echo "     -o /path/to/output"
}

# Function to log messages
log() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        "INFO")
            echo -e "${GREEN}[INFO]${NC} ${timestamp} - $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} ${timestamp} - $message"
            ;;
        "ERROR")
            echo -e "${RED}[ERROR]${NC} ${timestamp} - $message"
            ;;
        "DEBUG")
            if [ "$VERBOSE" = true ]; then
                echo -e "${BLUE}[DEBUG]${NC} ${timestamp} - $message"
            fi
            ;;
    esac
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            EQUIBIND_OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--protein)
            PROTEIN_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -e|--env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$EQUIBIND_OUTPUT_DIR" ] || [ -z "$PROTEIN_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    log "ERROR" "Missing required arguments"
    usage
    exit 1
fi

# Validate input files and directories
if [ ! -d "$EQUIBIND_OUTPUT_DIR" ]; then
    log "ERROR" "EquiBind output directory does not exist: $EQUIBIND_OUTPUT_DIR"
    exit 1
fi

if [ ! -f "$PROTEIN_FILE" ]; then
    log "ERROR" "Protein file does not exist: $PROTEIN_FILE"
    exit 1
fi

LIGANDS_DIR="$EQUIBIND_OUTPUT_DIR/ligands"
if [ ! -d "$LIGANDS_DIR" ]; then
    log "ERROR" "Ligands directory not found: $LIGANDS_DIR"
    exit 1
fi

SUCCESS_FILE="$LIGANDS_DIR/success.txt"
if [ ! -f "$SUCCESS_FILE" ]; then
    log "ERROR" "Success file not found: $SUCCESS_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
log "INFO" "Created output directory: $OUTPUT_DIR"

# Set paths
DELTA_LINF9_DIR="/vol/data/drug-design-pipeline/external/deltalinf9"
RESULTS_FILE="$OUTPUT_DIR/delta_linf9_results.csv"
RANKED_FILE="$OUTPUT_DIR/ligands_ranked.csv"
LOG_FILE="$OUTPUT_DIR/processing.log"

# Initialize results file
echo "ligand_id,ligand_file,xgb_score,linf9_score,processing_time" > "$RESULTS_FILE"

log "INFO" "Starting Delta LinF9 processing"
log "INFO" "EquiBind output directory: $EQUIBIND_OUTPUT_DIR"
log "INFO" "Protein file: $PROTEIN_FILE"
log "INFO" "Output directory: $OUTPUT_DIR"
log "INFO" "Results will be saved to: $RESULTS_FILE"

# Activate conda environment
log "INFO" "Activating conda environment: $CONDA_ENV"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Change to Delta LinF9 directory
cd "$DELTA_LINF9_DIR"

# Read successful ligand indices
successful_ligands=()
while IFS= read -r line; do
    if [[ "$line" =~ ^[0-9]+$ ]]; then
        successful_ligands+=("$line")
    fi
done < "$SUCCESS_FILE"

log "INFO" "Found ${#successful_ligands[@]} successful ligands to process"

# Process each ligand
total_ligands=${#successful_ligands[@]}
processed=0
failed=0

for ligand_id in "${successful_ligands[@]}"; do
    ligand_file="$LIGANDS_DIR/${ligand_id}.sdf"
    
    if [ ! -f "$ligand_file" ]; then
        log "WARN" "Ligand file not found: $ligand_file"
        ((failed++))
        continue
    fi
    
    log "DEBUG" "Processing ligand $ligand_id ($(($processed + 1))/$total_ligands)"
    
    # Record start time
    start_time=$(date +%s)
    
    # Run Delta LinF9
    temp_output=$(mktemp)
    if timeout 300 python script/runXGB.py "$PROTEIN_FILE" "$ligand_file" > "$temp_output" 2>&1; then
        # Extract XGB score
        xgb_score=$(grep "XGB (in pK)" "$temp_output" | sed 's/.*XGB (in pK) :\s*//' | tr -d ' ')
        
        # Try to extract LinF9 score (if available in output)
        linf9_score=$(grep -o "LinF9.*[0-9.-]\+" "$temp_output" | tail -1 | grep -o "[0-9.-]\+$" || echo "N/A")
        
        # Calculate processing time
        end_time=$(date +%s)
        processing_time=$((end_time - start_time))
        
        if [ -n "$xgb_score" ] && [ "$xgb_score" != "nan" ]; then
            echo "$ligand_id,$ligand_file,$xgb_score,$linf9_score,$processing_time" >> "$RESULTS_FILE"
            log "INFO" "Ligand $ligand_id: XGB score = $xgb_score"
            ((processed++))
        else
            log "WARN" "Failed to extract valid XGB score for ligand $ligand_id"
            echo "$ligand_id,$ligand_file,FAILED,$linf9_score,$processing_time" >> "$RESULTS_FILE"
            ((failed++))
        fi
    else
        log "ERROR" "Delta LinF9 processing failed for ligand $ligand_id"
        processing_time=$(($(date +%s) - start_time))
        echo "$ligand_id,$ligand_file,TIMEOUT,N/A,$processing_time" >> "$RESULTS_FILE"
        ((failed++))
    fi
    
    rm -f "$temp_output"
done

log "INFO" "Processing completed: $processed successful, $failed failed"

# Create ranked results (excluding failed entries)
log "INFO" "Creating ranked results..."

# Sort by XGB score in descending order (excluding FAILED and TIMEOUT entries)
{
    head -1 "$RESULTS_FILE"  # Header
    tail -n +2 "$RESULTS_FILE" | grep -v "FAILED\|TIMEOUT" | sort -t, -k3 -nr
    tail -n +2 "$RESULTS_FILE" | grep "FAILED\|TIMEOUT"  # Failed entries at the bottom
} > "$RANKED_FILE"

log "INFO" "Ranked results saved to: $RANKED_FILE"

# Generate summary report
SUMMARY_FILE="$OUTPUT_DIR/summary_report.txt"
{
    echo "Delta LinF9 Processing Summary Report"
    echo "====================================="
    echo "Generated: $(date)"
    echo ""
    echo "Input Parameters:"
    echo "  EquiBind Output: $EQUIBIND_OUTPUT_DIR"
    echo "  Protein File: $PROTEIN_FILE"
    echo "  Output Directory: $OUTPUT_DIR"
    echo ""
    echo "Processing Results:"
    echo "  Total ligands: $total_ligands"
    echo "  Successfully processed: $processed"
    echo "  Failed: $failed"
    echo ""
    echo "Top 10 Ligands by XGB Score:"
    echo "Rank | Ligand ID | XGB Score | LinF9 Score | Processing Time"
    echo "-----+-----------+-----------+-------------+----------------"
    
    # Show top 10 successful results
    tail -n +2 "$RANKED_FILE" | grep -v "FAILED\|TIMEOUT" | head -10 | nl -w4 -s" | " | \
    awk -F',' '{printf "%s %-9s | %-9s | %-11s | %s s\n", $1, $2, $3, $4, $5}'
    
} > "$SUMMARY_FILE"

log "INFO" "Summary report saved to: $SUMMARY_FILE"

# Display summary
echo ""
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}      Delta LinF9 Processing Complete      ${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
cat "$SUMMARY_FILE"
echo ""
echo -e "${BLUE}Files generated:${NC}"
echo "  📊 Results: $RESULTS_FILE"
echo "  🏆 Ranked: $RANKED_FILE"
echo "  📋 Summary: $SUMMARY_FILE"
echo ""

if [ $failed -gt 0 ]; then
    log "WARN" "Some ligands failed processing. Check the results files for details."
    exit 1
else
    log "INFO" "All ligands processed successfully!"
    exit 0
fi
