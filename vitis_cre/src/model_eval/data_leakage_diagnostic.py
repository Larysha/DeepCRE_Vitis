"""
Data Leakage Diagnostic Script for DeepCRE Cross-Validation
Identifies potential data leakage between standard and cross-prediction datasets.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from typing import Dict, Set, Tuple
import pickle


def load_validation_genes(pickle_path: str, pickle_key: str) -> Set[str]:
    """
    Load validation genes from pickle file.
    
    Args:
        pickle_path: Path to validation genes pickle file
        pickle_key: Key to access validation genes in pickle dict
    
    Returns:
        Set of validation gene IDs
    """
    try:
        with open(pickle_path, 'rb') as f:
            validation_data = pickle.load(f)
        
        if pickle_key in validation_data:
            return set(validation_data[pickle_key])
        else:
            print(f"Warning: Key '{pickle_key}' not found in pickle file.")
            print(f"Available keys: {list(validation_data.keys())}")
            return set()
    except FileNotFoundError:
        print(f"Warning: Pickle file not found at {pickle_path}")
        return set()


def analyze_gene_overlap(standard_file: str, cross_file: str, 
                        validation_pickle: str = None, pickle_key: str = None) -> Dict:
    """
    Comprehensive analysis of gene overlap between datasets.
    
    Args:
        standard_file: Path to standard predictions CSV
        cross_file: Path to cross-predictions CSV
        validation_pickle: Optional path to validation genes pickle
        pickle_key: Key for validation genes in pickle
    
    Returns:
        Dictionary containing overlap analysis results
    """
    print("Loading datasets...")
    
    # Load datasets
    std_df = pd.read_csv(standard_file)
    cross_df = pd.read_csv(cross_file)
    
    print(f"Standard dataset: {len(std_df)} genes")
    print(f"Cross dataset: {len(cross_df)} genes")
    
    # Get gene sets
    std_genes = set(std_df['genes'].tolist())
    cross_genes = set(cross_df['genes'].tolist())
    
    # Calculate overlaps
    overlap_genes = std_genes.intersection(cross_genes)
    std_only = std_genes - cross_genes
    cross_only = cross_genes - std_genes
    
    results = {
        'std_total': len(std_genes),
        'cross_total': len(cross_genes),
        'overlap_count': len(overlap_genes),
        'overlap_percent': len(overlap_genes) / len(std_genes) * 100,
        'std_only': len(std_only),
        'cross_only': len(cross_only),
        'overlap_genes': overlap_genes,
        'std_genes': std_genes,
        'cross_genes': cross_genes
    }
    
    # Load validation genes if provided
    if validation_pickle and pickle_key:
        val_genes = load_validation_genes(validation_pickle, pickle_key)
        if val_genes:
            val_in_std = val_genes.intersection(std_genes)
            val_in_cross = val_genes.intersection(cross_genes)
            val_in_both = val_genes.intersection(overlap_genes)
            
            results.update({
                'validation_total': len(val_genes),
                'val_in_std': len(val_in_std),
                'val_in_cross': len(val_in_cross),
                'val_in_both': len(val_in_both),
                'val_overlap_with_std_percent': len(val_in_std) / len(val_genes) * 100,
                'val_overlap_with_cross_percent': len(val_in_cross) / len(val_genes) * 100
            })
    
    return results


def analyze_class_distributions(standard_file: str, cross_file: str) -> Dict:
    """
    Compare class distributions between datasets.
    
    Args:
        standard_file: Path to standard predictions CSV
        cross_file: Path to cross-predictions CSV
    
    Returns:
        Dictionary with class distribution analysis
    """
    std_df = pd.read_csv(standard_file)
    cross_df = pd.read_csv(cross_file)
    
    # Filter to binary targets only
    std_binary = std_df[std_df['true_targets'].isin([0, 1])]
    cross_binary = cross_df[cross_df['true_targets'].isin([0, 1])]
    
    # Get class counts
    std_class_counts = std_binary['true_targets'].value_counts().sort_index()
    cross_class_counts = cross_binary['true_targets'].value_counts().sort_index()
    
    results = {
        'std_class_0': std_class_counts.get(0, 0),
        'std_class_1': std_class_counts.get(1, 0),
        'cross_class_0': cross_class_counts.get(0, 0),
        'cross_class_1': cross_class_counts.get(1, 0),
        'std_total_binary': len(std_binary),
        'cross_total_binary': len(cross_binary),
        'std_class_balance': std_class_counts.get(1, 0) / (std_class_counts.get(0, 0) + std_class_counts.get(1, 0)),
        'cross_class_balance': cross_class_counts.get(1, 0) / (cross_class_counts.get(0, 0) + cross_class_counts.get(1, 0))
    }
    
    return results


def analyze_chromosome_distribution(standard_file: str, cross_file: str) -> Dict:
    """
    Analyze gene distribution across chromosomes.
    
    Args:
        standard_file: Path to standard predictions CSV
        cross_file: Path to cross-predictions CSV
    
    Returns:
        Dictionary with chromosome distribution analysis
    """
    std_df = pd.read_csv(standard_file)
    cross_df = pd.read_csv(cross_file)
    
    # Extract chromosome from gene IDs
    def extract_chromosome(gene_id):
        # Assuming format like "Vitvi05_01chr01g00010"
        if 'chr' in gene_id:
            try:
                chrom_part = gene_id.split('chr')[1]
                chrom_num = chrom_part[:2]  # Get first 2 characters after 'chr'
                return f"chr{chrom_num}"
            except:
                return "unknown"
        return "unknown"
    
    std_df['chromosome'] = std_df['genes'].apply(extract_chromosome)
    cross_df['chromosome'] = cross_df['genes'].apply(extract_chromosome)
    
    std_chrom_counts = std_df['chromosome'].value_counts()
    cross_chrom_counts = cross_df['chromosome'].value_counts()
    
    return {
        'std_chrom_counts': std_chrom_counts,
        'cross_chrom_counts': cross_chrom_counts,
        'chromosomes_std': set(std_chrom_counts.index),
        'chromosomes_cross': set(cross_chrom_counts.index)
    }


def create_diagnostic_plots(overlap_results: Dict, class_results: Dict, 
                          chrom_results: Dict, output_dir: str = "diagnostic_results"):
    """
    Create comprehensive diagnostic visualizations.
    
    Args:
        overlap_results: Results from gene overlap analysis
        class_results: Results from class distribution analysis
        chrom_results: Results from chromosome analysis
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Data Leakage Diagnostic Analysis', fontsize=16, fontweight='bold')
    
    # 1. Gene Overlap Venn-style visualization
    ax1 = axes[0, 0]
    overlap_data = [
        overlap_results['std_only'],
        overlap_results['overlap_count'],
        overlap_results['cross_only']
    ]
    labels = ['Standard\nOnly', 'Overlap', 'Cross\nOnly']
    colors = ['#AF91AF', '#804a5f', '#478fca']
    
    bars = ax1.bar(labels, overlap_data, color=colors, alpha=0.8)
    ax1.set_title('Gene Set Overlap Analysis')
    ax1.set_ylabel('Number of Genes')
    
    # Add value labels
    for bar, value in zip(bars, overlap_data):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{value:,}\n({value/overlap_results["std_total"]*100:.1f}%)',
                ha='center', va='bottom', fontweight='bold')
    
    # 2. Class Distribution Comparison
    ax2 = axes[0, 1]
    class_data = np.array([
        [class_results['std_class_0'], class_results['std_class_1']],
        [class_results['cross_class_0'], class_results['cross_class_1']]
    ])
    
    x = np.arange(2)
    width = 0.35
    
    ax2.bar(x - width/2, class_data[0], width, label='Standard', color='#AF91AF', alpha=0.8)
    ax2.bar(x + width/2, class_data[1], width, label='Cross-Prediction', color='#804a5f', alpha=0.8)
    
    ax2.set_xlabel('Expression Class')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Class Distribution Comparison')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Low (0)', 'High (1)'])
    ax2.legend()
    
    # Add value labels
    for i, (std_val, cross_val) in enumerate(class_data.T):
        ax2.text(i - width/2, std_val + std_val*0.01, f'{std_val:,}', ha='center', va='bottom')
        ax2.text(i + width/2, cross_val + cross_val*0.01, f'{cross_val:,}', ha='center', va='bottom')
    
    # 3. Validation Gene Analysis (if available)
    ax3 = axes[0, 2]
    if 'validation_total' in overlap_results:
        val_data = [
            overlap_results['validation_total'] - overlap_results['val_in_std'] - overlap_results['val_in_cross'] + overlap_results['val_in_both'],
            overlap_results['val_in_std'] - overlap_results['val_in_both'],
            overlap_results['val_in_cross'] - overlap_results['val_in_both'],
            overlap_results['val_in_both']
        ]
        val_labels = ['Val Only', 'Val+Std', 'Val+Cross', 'Val+Both']
        ax3.pie(val_data, labels=val_labels, autopct='%1.1f%%', startangle=90)
        ax3.set_title('Validation Gene Distribution')
    else:
        ax3.text(0.5, 0.5, 'No validation\ngene data\nprovided', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Validation Gene Analysis')
    
    # 4. Chromosome Distribution Standard
    ax4 = axes[1, 0]
    std_chroms = chrom_results['std_chrom_counts']
    if len(std_chroms) > 0:
        ax4.bar(range(len(std_chroms)), std_chroms.values, color='#AF91AF', alpha=0.8)
        ax4.set_xlabel('Chromosome')
        ax4.set_ylabel('Gene Count')
        ax4.set_title('Standard Dataset: Genes per Chromosome')
        ax4.set_xticks(range(len(std_chroms)))
        ax4.set_xticklabels(std_chroms.index, rotation=45)
    
    # 5. Chromosome Distribution Cross
    ax5 = axes[1, 1]
    cross_chroms = chrom_results['cross_chrom_counts']
    if len(cross_chroms) > 0:
        ax5.bar(range(len(cross_chroms)), cross_chroms.values, color='#804a5f', alpha=0.8)
        ax5.set_xlabel('Chromosome')
        ax5.set_ylabel('Gene Count')
        ax5.set_title('Cross-Prediction Dataset: Genes per Chromosome')
        ax5.set_xticks(range(len(cross_chroms)))
        ax5.set_xticklabels(cross_chroms.index, rotation=45)
    
    # 6. Summary Statistics
    ax6 = axes[1, 2]
    ax6.axis('off')
    
    summary_text = f"""
    DATA LEAKAGE DIAGNOSTIC SUMMARY
    
    Gene Overlap:
    • Total overlap: {overlap_results['overlap_count']:,} genes
    • Overlap percentage: {overlap_results['overlap_percent']:.1f}%
    • Standard only: {overlap_results['std_only']:,}
    • Cross only: {overlap_results['cross_only']:,}
    
    Class Balance Similarity:
    • Standard: {class_results['std_class_balance']:.3f}
    • Cross: {class_results['cross_class_balance']:.3f}
    • Difference: {abs(class_results['std_class_balance'] - class_results['cross_class_balance']):.3f}
    
    LEAKAGE RISK ASSESSMENT:
    """
    
    # Risk assessment
    if overlap_results['overlap_percent'] > 90:
        risk_text = "🚨 CRITICAL: >90% gene overlap"
    elif overlap_results['overlap_percent'] > 50:
        risk_text = "⚠️  HIGH: >50% gene overlap"
    elif overlap_results['overlap_percent'] > 10:
        risk_text = "⚠️  MODERATE: >10% gene overlap"
    else:
        risk_text = "✅ LOW: <10% gene overlap"
    
    summary_text += risk_text
    
    ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white", 
                     edgecolor='#AF91AF', linewidth=1.5, alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/data_leakage_diagnostic.png", dpi=300, bbox_inches='tight')
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Diagnose data leakage in DeepCRE cross-validation')
    parser.add_argument('--standard', '-s', required=True,
                       help='Path to standard predictions CSV file')
    parser.add_argument('--cross', '-c', required=True,
                       help='Path to cross-predictions CSV file')
    parser.add_argument('--validation_pickle', '-v', required=False,
                       help='Path to validation genes pickle file')
    parser.add_argument('--pickle_key', '-k', required=False, default='vitis',
                       help='Key for validation genes in pickle (default: vitis)')
    parser.add_argument('--output', '-o', default='diagnostic_results',
                       help='Output directory for diagnostic results')
    
    args = parser.parse_args()
    
    print("="*60)
    print("DATA LEAKAGE DIAGNOSTIC ANALYSIS")
    print("="*60)
    
    # Run analyses
    overlap_results = analyze_gene_overlap(args.standard, args.cross, 
                                          args.validation_pickle, args.pickle_key)
    class_results = analyze_class_distributions(args.standard, args.cross)
    chrom_results = analyze_chromosome_distribution(args.standard, args.cross)
    
    # Print results
    print(f"\nGENE OVERLAP ANALYSIS:")
    print(f"Standard dataset: {overlap_results['std_total']:,} genes")
    print(f"Cross dataset: {overlap_results['cross_total']:,} genes")
    print(f"Overlap: {overlap_results['overlap_count']:,} genes ({overlap_results['overlap_percent']:.1f}%)")
    print(f"Standard only: {overlap_results['std_only']:,} genes")
    print(f"Cross only: {overlap_results['cross_only']:,} genes")
    
    print(f"\nCLASS DISTRIBUTION ANALYSIS:")
    print(f"Standard - Class 0: {class_results['std_class_0']:,}, Class 1: {class_results['std_class_1']:,}")
    print(f"Cross - Class 0: {class_results['cross_class_0']:,}, Class 1: {class_results['cross_class_1']:,}")
    print(f"Class balance difference: {abs(class_results['std_class_balance'] - class_results['cross_class_balance']):.4f}")
    
    if 'validation_total' in overlap_results:
        print(f"\nVALIDATION GENE ANALYSIS:")
        print(f"Total validation genes: {overlap_results['validation_total']:,}")
        print(f"Validation genes in standard: {overlap_results['val_in_std']:,} ({overlap_results['val_overlap_with_std_percent']:.1f}%)")
        print(f"Validation genes in cross: {overlap_results['val_in_cross']:,} ({overlap_results['val_overlap_with_cross_percent']:.1f}%)")
        print(f"Validation genes in both: {overlap_results['val_in_both']:,}")
    
    # Risk assessment
    print(f"\nRISK ASSESSMENT:")
    if overlap_results['overlap_percent'] > 90:
        print("🚨 CRITICAL DATA LEAKAGE: >90% gene overlap detected!")
        print("   Cross-validation results are invalid.")
    elif overlap_results['overlap_percent'] > 50:
        print("⚠️  HIGH LEAKAGE RISK: >50% gene overlap detected.")
        print("   Results should be interpreted with extreme caution.")
    elif overlap_results['overlap_percent'] > 10:
        print("⚠️  MODERATE LEAKAGE RISK: >10% gene overlap detected.")
        print("   Some performance inflation expected.")
    else:
        print("✅ LOW LEAKAGE RISK: <10% gene overlap.")
        print("   Cross-validation appears properly isolated.")
    
    # Check class balance similarity
    balance_diff = abs(class_results['std_class_balance'] - class_results['cross_class_balance'])
    if balance_diff < 0.001:
        print("⚠️  SUSPICIOUS: Identical class balances suggest same dataset.")
    
    # Create diagnostic plots
    create_diagnostic_plots(overlap_results, class_results, chrom_results, args.output)
    
    # Save detailed results
    output_path = Path(args.output)
    output_path.mkdir(exist_ok=True)
    
    # Save overlap gene lists for further investigation
    with open(output_path / "overlapping_genes.txt", 'w') as f:
        f.write("Genes present in both standard and cross-prediction datasets:\n")
        f.write("="*50 + "\n")
        for gene in sorted(overlap_results['overlap_genes']):
            f.write(f"{gene}\n")
    
    print(f"\nDiagnostic results saved to {args.output}/")
    print("Files created:")
    print("- data_leakage_diagnostic.png")
    print("- overlapping_genes.txt")


if __name__ == "__main__":
    main()