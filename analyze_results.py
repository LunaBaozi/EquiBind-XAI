#!/usr/bin/env python3
"""
Delta LinF9 Results Analyzer
Helper script for analyzing and visualizing Delta LinF9 XGB scoring results
"""

import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def load_results(results_file):
    """Load and validate the results CSV file"""
    try:
        df = pd.read_csv(results_file)
        # Convert XGB scores to numeric, handling FAILED/TIMEOUT entries
        df['xgb_score_numeric'] = pd.to_numeric(df['xgb_score'], errors='coerce')
        return df
    except Exception as e:
        print(f"Error loading results file: {e}")
        sys.exit(1)

def generate_summary_stats(df):
    """Generate summary statistics for the XGB scores"""
    valid_scores = df['xgb_score_numeric'].dropna()
    
    if len(valid_scores) == 0:
        return "No valid XGB scores found."
    
    stats = {
        'count': len(valid_scores),
        'mean': valid_scores.mean(),
        'std': valid_scores.std(),
        'min': valid_scores.min(),
        'max': valid_scores.max(),
        'median': valid_scores.median(),
        'q25': valid_scores.quantile(0.25),
        'q75': valid_scores.quantile(0.75)
    }
    
    return stats

def create_visualizations(df, output_dir):
    """Create visualizations of the XGB scores"""
    valid_df = df.dropna(subset=['xgb_score_numeric'])
    
    if len(valid_df) == 0:
        print("No valid scores for visualization")
        # Create empty plot file
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, 'No valid XGB scores to plot', ha='center', va='center', fontsize=14)
        plt.title('Delta LinF9 XGB Score Analysis - No Data')
        plot_file = os.path.join(output_dir, 'xgb_analysis.png')
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Empty plot saved to: {plot_file}")
        return
    
    # Set up the plotting style - use a simpler style
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Delta LinF9 XGB Score Analysis', fontsize=16, fontweight='bold')
    
    # 1. Histogram of XGB scores
    axes[0, 0].hist(valid_df['xgb_score_numeric'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].set_xlabel('XGB Score (pK)')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of XGB Scores')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Box plot
    axes[0, 1].boxplot(valid_df['xgb_score_numeric'])
    axes[0, 1].set_ylabel('XGB Score (pK)')
    axes[0, 1].set_title('XGB Score Box Plot')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Scatter plot: Ligand ID vs XGB Score
    axes[1, 0].scatter(valid_df['ligand_id'], valid_df['xgb_score_numeric'], alpha=0.6, color='coral')
    axes[1, 0].set_xlabel('Ligand ID')
    axes[1, 0].set_ylabel('XGB Score (pK)')
    axes[1, 0].set_title('XGB Score by Ligand ID')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Processing time vs XGB score
    if 'processing_time' in valid_df.columns:
        valid_time_df = valid_df.dropna(subset=['processing_time'])
        if len(valid_time_df) > 0:
            axes[1, 1].scatter(valid_time_df['processing_time'], valid_time_df['xgb_score_numeric'], 
                             alpha=0.6, color='lightgreen')
            axes[1, 1].set_xlabel('Processing Time (seconds)')
            axes[1, 1].set_ylabel('XGB Score (pK)')
            axes[1, 1].set_title('XGB Score vs Processing Time')
            axes[1, 1].grid(True, alpha=0.3)
        else:
            axes[1, 1].text(0.5, 0.5, 'No valid processing time data', 
                           ha='center', va='center', transform=axes[1, 1].transAxes)
    else:
        axes[1, 1].text(0.5, 0.5, 'No processing time data available', 
                       ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    
    # Save the plot
    plot_file = os.path.join(output_dir, 'xgb_analysis.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Visualization saved to: {plot_file}")

def export_top_ligands(df, output_dir, top_n=10):
    """Export the top N ligands to a separate file"""
    valid_df = df.dropna(subset=['xgb_score_numeric'])
    
    output_file = os.path.join(output_dir, f'top_{top_n}_ligands.csv')
    
    if len(valid_df) > 0:
        top_ligands = valid_df.nlargest(top_n, 'xgb_score_numeric')
        top_ligands.to_csv(output_file, index=False)
        print(f"Top {top_n} ligands exported to: {output_file}")
        return top_ligands
    else:
        # Create empty file with headers
        empty_df = pd.DataFrame(columns=['ligand_id', 'xgb_score', 'processing_time', 'status', 'ligand_file'])
        empty_df.to_csv(output_file, index=False)
        print(f"No valid ligands found. Empty file created: {output_file}")
        return empty_df

def main():
    parser = argparse.ArgumentParser(description='Analyze Delta LinF9 XGB scoring results')
    parser.add_argument('results_file', help='Path to the results CSV file')
    parser.add_argument('-o', '--output', default='.', help='Output directory for analysis files')
    parser.add_argument('-t', '--top', type=int, default=10, help='Number of top ligands to export')
    parser.add_argument('--plot', action='store_true', help='Generate visualization plots')
    parser.add_argument('--stats', action='store_true', help='Show detailed statistics')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.results_file):
        print(f"Error: Results file not found: {args.results_file}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Load results
    print(f"Loading results from: {args.results_file}")
    df = load_results(args.results_file)
    
    print(f"Total ligands: {len(df)}")
    print(f"Successfully scored: {len(df.dropna(subset=['xgb_score_numeric']))}")
    print(f"Failed: {len(df) - len(df.dropna(subset=['xgb_score_numeric']))}")
    
    # Generate statistics
    if args.stats:
        print("\n" + "="*50)
        print("STATISTICAL SUMMARY")
        print("="*50)
        stats = generate_summary_stats(df)
        if isinstance(stats, dict):
            for key, value in stats.items():
                print(f"{key.upper():12}: {value:.4f}")
        else:
            print(stats)
    
    # Export top ligands
    print(f"\nExporting top {args.top} ligands...")
    top_ligands = export_top_ligands(df, args.output, args.top)
    
    if len(top_ligands) > 0:
        print(f"\nTOP {args.top} LIGANDS:")
        print("-" * 60)
        for idx, row in top_ligands.iterrows():
            print(f"{row['ligand_id']:3} | {row['xgb_score_numeric']:8.3f} | {row.get('ligand_file', 'N/A')}")
    else:
        print("No valid ligands found for ranking.")
    
    # Generate plots (always create the PNG file)
    print(f"\nGenerating visualizations...")
    try:
        create_visualizations(df, args.output)
    except ImportError:
        print("Warning: matplotlib/seaborn not available. Creating empty plot file.")
        # Create empty plot file
        plot_file = os.path.join(args.output, 'xgb_analysis.png')
        with open(plot_file, 'w') as f:
            f.write("# Empty plot file - matplotlib not available\n")
    except Exception as e:
        print(f"Error generating visualizations: {e}")
        # Create empty plot file
        plot_file = os.path.join(args.output, 'xgb_analysis.png')
        with open(plot_file, 'w') as f:
            f.write("# Empty plot file - error occurred\n")
    
    print(f"\nAnalysis complete. Output files saved to: {args.output}")

if __name__ == "__main__":
    main()
