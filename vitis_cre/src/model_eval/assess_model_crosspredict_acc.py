"""
DeepCRE Chromosomal Cross-Prediction Analysis
Analyzes performance of chromosome-specific models on all chromosomes.
Creates heatmaps and performance metrics
"""

import pandas as pd
import numpy as np
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score, 
    roc_auc_score, average_precision_score, confusion_matrix,
    classification_report
)
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re
from typing import Dict, List, Tuple, Optional


# Custom colour palette matching your theme
CUSTOM_COLORS = ["#e1c7cd", "#671436", "#0C439B", "#1d778b"]
EXTENDED_COLORS = CUSTOM_COLORS + ['#a8dadc', '#457b9d', '#1d3557', '#f1faee']


def parse_model_columns(df: pd.DataFrame) -> List[str]:
    """
    Parse model column names to extract chromosome information.
    
    Args:
        df: DataFrame with cross-prediction results
        
    Returns:
        List of chromosome model column names
    """
    model_columns = []
    for col in df.columns:
        # Look for model files ending in .h5
        if col.endswith('.h5') and 'chr' in col:
            model_columns.append(col)
    
    # Sort chromosomes numerically
    def extract_chr_num(col_name):
        match = re.search(r'chr(\d+)', col_name)
        return int(match.group(1)) if match else 0
    
    model_columns.sort(key=extract_chr_num)
    return model_columns


def extract_chromosome_from_gene(gene_id: str) -> Optional[str]:
    """
    Extract chromosome information from gene ID.
    
    Args:
        gene_id: Gene identifier (e.g., 'Vitvi05_01chr01g00010')
        
    Returns:
        Chromosome identifier or None if not found
    """
    # Look for chromosome pattern in gene ID
    match = re.search(r'chr(\d+)', gene_id)
    if match:
        return f"chr{match.group(1).zfill(2)}"  # Ensure 2-digit format
    return None


def calculate_cross_prediction_metrics(df: pd.DataFrame, model_columns: List[str], 
                                     threshold: float = 0.5) -> pd.DataFrame:
    """
    Calculate cross-prediction metrics for each model-chromosome combination.
    
    Args:
        df: DataFrame with cross-prediction results
        model_columns: List of model column names
        threshold: Prediction threshold for binary classification
        
    Returns:
        DataFrame with metrics for each model-chromosome pair
    """
    results = []
    
    # Add chromosome column to dataframe
    df['gene_chromosome'] = df['genes'].apply(extract_chromosome_from_gene)
    
    # Get unique chromosomes from genes
    chromosomes = sorted(df['gene_chromosome'].dropna().unique())
    
    for model_col in model_columns:
        # Extract training chromosome from model name
        train_chr_match = re.search(r'chr(\d+)', model_col)
        train_chr = f"chr{train_chr_match.group(1).zfill(2)}" if train_chr_match else "unknown"
        
        for test_chr in chromosomes:
            # Filter genes from this test chromosome
            chr_genes = df[df['gene_chromosome'] == test_chr].copy()
            
            if len(chr_genes) == 0:
                continue
                
            # Get predictions for this model
            y_true = chr_genes['true_targets']
            y_pred_probs = chr_genes[model_col]
            y_pred_binary = (y_pred_probs >= threshold).astype(int)
            
            # Calculate metrics
            accuracy = accuracy_score(y_true, y_pred_binary)
            precision = precision_score(y_true, y_pred_binary, zero_division=0)
            recall = recall_score(y_true, y_pred_binary, zero_division=0)
            f1 = f1_score(y_true, y_pred_binary, zero_division=0)
            
            # Only calculate AUC if we have both classes
            if len(np.unique(y_true)) == 2:
                roc_auc = roc_auc_score(y_true, y_pred_probs)
                pr_auc = average_precision_score(y_true, y_pred_probs)
            else:
                roc_auc = np.nan
                pr_auc = np.nan
            
            results.append({
                'train_chromosome': train_chr,
                'test_chromosome': test_chr,
                'model_column': model_col,
                'n_genes': len(chr_genes),
                'accuracy': accuracy,
                'precision': precision,
                'recall': recall,
                'f1_score': f1,
                'roc_auc': roc_auc,
                'pr_auc': pr_auc,
                'is_self_prediction': train_chr == test_chr
            })
    
    return pd.DataFrame(results)


def calculate_ensemble_metrics(df: pd.DataFrame, model_columns: List[str], 
                             threshold: float = 0.5) -> Dict:
    """
    Calculate ensemble prediction metrics using all models.
    
    Args:
        df: DataFrame with cross-prediction results
        model_columns: List of model column names
        threshold: Prediction threshold for binary classification
        
    Returns:
        Dictionary with ensemble metrics
    """
    # Calculate ensemble predictions (mean across all models)
    df['ensemble_pred_probs'] = df[model_columns].mean(axis=1)
    df['ensemble_pred_binary'] = (df['ensemble_pred_probs'] >= threshold).astype(int)
    
    # Calculate ensemble metrics
    y_true = df['true_targets']
    y_pred_probs = df['ensemble_pred_probs']
    y_pred_binary = df['ensemble_pred_binary']
    
    metrics = {
        'accuracy': accuracy_score(y_true, y_pred_binary),
        'accuracy_percent': accuracy_score(y_true, y_pred_binary) * 100,
        'precision': precision_score(y_true, y_pred_binary, zero_division=0),
        'recall': recall_score(y_true, y_pred_binary, zero_division=0),
        'f1_score': f1_score(y_true, y_pred_binary, zero_division=0),
        'n_genes': len(df)
    }
    
    if len(np.unique(y_true)) == 2:
        metrics['roc_auc'] = roc_auc_score(y_true, y_pred_probs)
        metrics['pr_auc'] = average_precision_score(y_true, y_pred_probs)
    
    return metrics


def create_heatmap(metrics_df: pd.DataFrame, metric: str = 'accuracy', 
                   output_dir: str = "vitis_cre/src/results/model_eval/cross_predict") -> None:
    """
    Create chromosomal cross-prediction heatmap.
    
    Args:
        metrics_df: DataFrame with cross-prediction metrics
        metric: Metric to visualise in heatmap
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Create pivot table for heatmap
    heatmap_data = metrics_df.pivot(index='train_chromosome', 
                                   columns='test_chromosome', 
                                   values=metric)
    
    # Sort chromosomes numerically
    def sort_chr(chr_list):
        return sorted(chr_list, key=lambda x: int(re.search(r'(\d+)', x).group(1)))
    
    sorted_chromosomes = sort_chr(heatmap_data.index.tolist())
    heatmap_data = heatmap_data.reindex(index=sorted_chromosomes, columns=sorted_chromosomes)
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Use a blue-white colour palette
    cmap = plt.cm.Blues
    
    # Add annotations
    im = ax.imshow(heatmap_data.values, cmap=cmap, aspect='auto')
    
    # Add colourbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(f'{metric.replace("_", " ").title()}', rotation=270, labelpad=20)
    
    # Set ticks and labels
    ax.set_xticks(range(len(heatmap_data.columns)))
    ax.set_yticks(range(len(heatmap_data.index)))
    ax.set_xticklabels([col.replace('chr', '') for col in heatmap_data.columns], rotation=45)
    ax.set_yticklabels([idx.replace('chr', '') for idx in heatmap_data.index])
    
    # Add text annotations
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            value = heatmap_data.iloc[i, j]
            if not np.isnan(value):
                # Highlight diagonal (self-predictions) with different colour
                text_color = 'white' if i == j else 'white' if value < 0.5 else 'black'
                ax.text(j, i, f'{value:.2f}', ha='center', va='center', 
                       color=text_color, fontweight='bold' if i == j else 'normal',
                       fontsize=9)
    
    # Labels and title
    ax.set_xlabel('Test Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('Training Chromosome', fontsize=12, fontweight='bold')
    ax.set_title(f'DeepCRE Chromosomal Cross-Prediction Performance\n{metric.replace("_", " ").title()} Matrix', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add subtitle explaining methodology
    fig.text(0.5, 0.95, 'Each model trained on one chromosome, tested on all chromosomes\nDiagonal shows self-prediction performance', 
             ha='center', va='top', fontsize=10, style='italic')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cross_prediction_heatmap_{metric}.png", dpi=300, bbox_inches='tight')
    plt.show()


def create_ensemble_evaluation_plots(df: pd.DataFrame, ensemble_metrics: Dict, model_columns: List[str],
                                   output_dir: str = "vitis_cre/src/results/model_eval/cross_predict") -> None:
    """
    Create ensemble-focused evaluation plots similar to self-prediction analysis.
    
    Args:
        df: DataFrame with cross-prediction results
        ensemble_metrics: Dictionary with ensemble metrics
        model_columns: List of model column names
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Add main title explaining ensemble methodology
    fig.suptitle('DeepCRE Ensemble Performance Assessment\nAverage predictions across all 18 chromosome-specific models', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 1. Confusion Matrix
    ax1 = axes[0, 0]
    y_true = df['true_targets']
    y_pred = df['ensemble_pred_binary']
    
    cm = confusion_matrix(y_true, y_pred)
    cmap = sns.blend_palette([CUSTOM_COLORS[0], CUSTOM_COLORS[1]], as_cmap=True)
    sns.heatmap(cm, annot=True, fmt='d', cmap=cmap, ax=ax1, 
                xticklabels=['Predicted 0', 'Predicted 1'],
                yticklabels=['True 0', 'True 1'])
    ax1.set_title('Ensemble Confusion Matrix')
    ax1.set_ylabel('True Labels')
    ax1.set_xlabel('Predicted Labels')
    
    # 2. Probability Distribution by Class
    ax2 = axes[0, 1]
    for i, class_label in enumerate([0, 1]):
        class_probs = df[df['true_targets'] == class_label]['ensemble_pred_probs']
        ax2.hist(class_probs, bins=30, alpha=0.5, 
                label=f'True Class {class_label}', density=True, 
                color=CUSTOM_COLORS[i])
    ax2.axvline(x=0.5, color=CUSTOM_COLORS[3], linestyle='--', linewidth=2, label='Threshold (0.5)')
    ax2.set_xlabel('Ensemble Predicted Probability')
    ax2.set_ylabel('Density')
    ax2.set_title('Ensemble Probability Distribution by True Class')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Class Balance
    ax3 = axes[1, 0]
    class_counts = df['true_targets'].value_counts().sort_index()
    bars = ax3.bar(['Class 0\n(Low Expression)', 'Class 1\n(High Expression)'], 
                   class_counts.values, color=[CUSTOM_COLORS[0], CUSTOM_COLORS[1]])
    ax3.set_title('Class Distribution in Dataset')
    ax3.set_ylabel('Number of Genes')
    
    # Add value labels on bars
    for bar, count in zip(bars, class_counts.values):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                f'{count:,}', ha='center', va='bottom', fontweight='bold')
    
    # 4. Ensemble Metrics Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Calculate confusion matrix components for ensemble
    if cm.shape == (2, 2):
        tn, fp, fn, tp = cm.ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    else:
        tn = fp = fn = tp = 0
        specificity = sensitivity = 0
    
    # Create ensemble metrics text
    metrics_text = f"""
    ENSEMBLE PERFORMANCE METRICS
    
    Overall Accuracy: {ensemble_metrics['accuracy']:.3f} ({ensemble_metrics['accuracy_percent']:.1f}%)
    Precision: {ensemble_metrics['precision']:.3f}
    Recall (Sensitivity): {ensemble_metrics['recall']:.3f}
    Specificity: {specificity:.3f}
    F1-Score: {ensemble_metrics['f1_score']:.3f}
    
    ROC-AUC: {ensemble_metrics.get('roc_auc', 0):.3f}
    PR-AUC: {ensemble_metrics.get('pr_auc', 0):.3f}
    
    True Positives: {tp:,}
    True Negatives: {tn:,}
    False Positives: {fp:,}
    False Negatives: {fn:,}
    
    Total Genes: {ensemble_metrics['n_genes']:,}
    """
    
    # White background with subtle border
    ax4.text(0.1, 0.9, metrics_text, transform=ax4.transAxes, fontsize=12,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white", 
                      edgecolor=CUSTOM_COLORS[0], linewidth=1.5, alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/ensemble_performance_summary.png", dpi=300, bbox_inches='tight')
    plt.show()


def ensemble_threshold_analysis(df: pd.DataFrame, output_dir: str = "vitis_cre/src/results/model_eval/cross_predict") -> pd.DataFrame:
    """
    Analyse ensemble performance across different probability thresholds.
    
    Args:
        df: DataFrame with ensemble predictions
        output_dir: Directory to save plots
        
    Returns:
        DataFrame with threshold analysis results
    """
    thresholds = np.arange(0.1, 0.95, 0.05)
    threshold_metrics = []
    
    for threshold in thresholds:
        y_pred_thresh = (df['ensemble_pred_probs'] >= threshold).astype(int)
        
        acc = accuracy_score(df['true_targets'], y_pred_thresh)
        prec = precision_score(df['true_targets'], y_pred_thresh, zero_division=0)
        rec = recall_score(df['true_targets'], y_pred_thresh, zero_division=0)
        f1 = f1_score(df['true_targets'], y_pred_thresh, zero_division=0)
        
        threshold_metrics.append({
            'threshold': threshold,
            'accuracy': acc,
            'precision': prec,
            'recall': rec,
            'f1_score': f1
        })
    
    threshold_df = pd.DataFrame(threshold_metrics)
    
    # Plot threshold analysis with custom colours
    plt.figure(figsize=(12, 8))
    
    # Add figure title explaining ensemble methodology
    plt.suptitle('DeepCRE Ensemble Threshold Analysis\nPerformance across probability thresholds (averaged predictions)', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.subplot(2, 1, 1)
    plt.plot(threshold_df['threshold'], threshold_df['accuracy'], 'o-', 
             label='Accuracy', linewidth=2, color=CUSTOM_COLORS[0], markersize=6)
    plt.plot(threshold_df['threshold'], threshold_df['f1_score'], 's-', 
             label='F1-Score', linewidth=2, color=CUSTOM_COLORS[1], markersize=6)
    plt.axvline(x=0.5, color=CUSTOM_COLORS[3], linestyle='--', alpha=0.8, 
                linewidth=2, label='Default Threshold')
    plt.xlabel('Probability Threshold')
    plt.ylabel('Score')
    plt.title('Ensemble Performance vs Probability Threshold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 1, 2)
    plt.plot(threshold_df['threshold'], threshold_df['precision'], '^-', 
             label='Precision', linewidth=2, color=CUSTOM_COLORS[2], markersize=6)
    plt.plot(threshold_df['threshold'], threshold_df['recall'], 'v-', 
             label='Recall', linewidth=2, color=CUSTOM_COLORS[0], markersize=6)
    plt.axvline(x=0.5, color=CUSTOM_COLORS[3], linestyle='--', alpha=0.8, 
                linewidth=2, label='Default Threshold')
    plt.xlabel('Probability Threshold')
    plt.ylabel('Score')
    plt.title('Precision vs Recall Trade-off')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/ensemble_threshold_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Find optimal threshold (max F1-score)
    optimal_idx = threshold_df['f1_score'].idxmax()
    optimal_threshold = threshold_df.loc[optimal_idx, 'threshold']
    optimal_f1 = threshold_df.loc[optimal_idx, 'f1_score']
    
    print(f"\nOptimal threshold for ensemble F1-score: {optimal_threshold:.2f} (F1 = {optimal_f1:.3f})")
    
    return threshold_df


def main():
    parser = argparse.ArgumentParser(description='Analyse DeepCRE chromosomal cross-predictions')
    parser.add_argument('--predictions', '-p', required=True,
                       help='Path to cross-predictions CSV file')
    parser.add_argument('--threshold', '-t', type=float, default=0.5,
                       help='Probability threshold for binary classification (default: 0.5)')
    parser.add_argument('--output', '-o', default='vitis_cre/src/results/model_eval/cross_predict',
                       help='Output directory for plots and results')
    parser.add_argument('--metric', '-m', default='accuracy',
                       choices=['accuracy', 'precision', 'recall', 'f1_score', 'roc_auc', 'pr_auc'],
                       help='Metric to display in heatmap (default: accuracy)')
    
    args = parser.parse_args()
    
    print("Loading cross-prediction data...")
    df = pd.read_csv(args.predictions)
    
    # Parse model columns
    model_columns = parse_model_columns(df)
    print(f"Found {len(model_columns)} chromosome models")
    
    # Filter out medium expression genes (target = 2)
    df_binary = df[df['true_targets'].isin([0, 1])].copy()
    print(f"Using {len(df_binary)} genes with binary expression labels (excluded medium expression)")
    
    # Calculate cross-prediction metrics
    print("Calculating cross-prediction metrics...")
    metrics_df = calculate_cross_prediction_metrics(df_binary, model_columns, args.threshold)
    
    # Calculate ensemble metrics
    print("Calculating ensemble metrics...")
    ensemble_metrics = calculate_ensemble_metrics(df_binary, model_columns, args.threshold)
    
    # Print summary statistics
    print("\n" + "="*60)
    print("CHROMOSOMAL CROSS-PREDICTION SUMMARY")
    print("="*60)
    
    self_pred = metrics_df[metrics_df['is_self_prediction'] == True]
    cross_pred = metrics_df[metrics_df['is_self_prediction'] == False]
    
    print(f"Self-Prediction Performance (diagonal):")
    print(f"  Mean Accuracy: {self_pred['accuracy'].mean():.3f} ± {self_pred['accuracy'].std():.3f}")
    print(f"  Range: {self_pred['accuracy'].min():.3f} - {self_pred['accuracy'].max():.3f}")
    
    print(f"\nCross-Prediction Performance (off-diagonal):")
    print(f"  Mean Accuracy: {cross_pred['accuracy'].mean():.3f} ± {cross_pred['accuracy'].std():.3f}")
    print(f"  Range: {cross_pred['accuracy'].min():.3f} - {cross_pred['accuracy'].max():.3f}")
    
    print(f"\nEnsemble Performance:")
    print(f"  Accuracy: {ensemble_metrics['accuracy']:.3f} ({ensemble_metrics['accuracy_percent']:.1f}%)")
    print(f"  F1-Score: {ensemble_metrics['f1_score']:.3f}")
    print(f"  ROC-AUC: {ensemble_metrics.get('roc_auc', 'N/A'):.3f}")
    
    print("="*60)
    
    # Create visualisations
    print("\nGenerating heatmap...")
    create_heatmap(metrics_df, args.metric, args.output)
    
    print("\nGenerating ensemble evaluation plots...")
    create_ensemble_evaluation_plots(df_binary, ensemble_metrics, model_columns, args.output)
    
    print("\nGenerating threshold analysis...")
    threshold_df = ensemble_threshold_analysis(df_binary, args.output)
    
    # Save detailed results
    print("\nSaving detailed results...")
    metrics_df.to_csv(f"{args.output}/cross_prediction_metrics.csv", index=False)
    
    # Save ensemble predictions and threshold analysis
    df_binary[['genes', 'true_targets', 'ensemble_pred_probs', 'ensemble_pred_binary']].to_csv(
        f"{args.output}/ensemble_predictions.csv", index=False)
    
    threshold_df.to_csv(f"{args.output}/ensemble_threshold_analysis.csv", index=False)
    
    # Save ensemble classification report
    report = classification_report(df_binary['true_targets'], df_binary['ensemble_pred_binary'], 
                                 target_names=['Low Expression', 'High Expression'])
    
    with open(f"{args.output}/ensemble_classification_report.txt", 'w') as f:
        f.write("DeepCRE Ensemble Model Classification Report\n")
        f.write("Average predictions across 18 chromosome-specific models\n")
        f.write("="*60 + "\n\n")
        f.write(report)
    
    print(f"\nResults saved to {args.output}/")
    print("Files created:")
    print("- cross_prediction_heatmap_[metric].png (diagnostic)")
    print("- ensemble_performance_summary.png")
    print("- ensemble_threshold_analysis.png")
    print("- ensemble_predictions.csv")
    print("- ensemble_classification_report.txt")
    print("- cross_prediction_metrics.csv")
    print("- ensemble_threshold_analysis.csv")


if __name__ == "__main__":
    main()