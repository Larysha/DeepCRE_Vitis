"""
Model accuracy assessment for DeepCRE self-predictions
Calculates comprehensive accuracy metrics for binary classification model outputs.
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


# Custom color palette
CUSTOM_COLORS = ["#a085a1", '#5b2339', "#5777ac", '#8d77ab' ]
EXTENDED_COLORS = CUSTOM_COLORS + ['#a8dadc', '#457b9d', '#1d3557', '#f1faee']


def load_predictions(file_path: str, threshold: float = 0.5) -> pd.DataFrame:
    """
    Load predictions and convert probabilities to binary predictions.
    
    Args:
        file_path: Path to CSV file with columns: true_targets, pred_probs, genes
        threshold: Probability threshold for binary classification (default: 0.5)
    
    Returns:
        DataFrame with original data plus binary predictions
    """
    df = pd.read_csv(file_path)
    
    # Validate required columns
    required_cols = ['true_targets', 'pred_probs', 'genes']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    
    # Convert probabilities to binary predictions
    df['pred_binary'] = (df['pred_probs'] >= threshold).astype(int)
    
    # Filter out any non-binary targets (if present)
    df = df[df['true_targets'].isin([0, 1])]
    
    print(f"Loaded {len(df)} predictions")
    print(f"Class distribution - 0: {sum(df['true_targets'] == 0)}, 1: {sum(df['true_targets'] == 1)}")
    
    return df


def calculate_metrics(y_true, y_pred, y_probs) -> dict:
    """
    Calculate comprehensive classification metrics.
    
    Args:
        y_true: True binary labels
        y_pred: Predicted binary labels  
        y_probs: Predicted probabilities
    
    Returns:
        Dictionary of metric values
    """
    metrics = {}
    
    # Basic classification metrics
    metrics['accuracy'] = accuracy_score(y_true, y_pred)
    metrics['accuracy_percent'] = metrics['accuracy'] * 100  # Add percentage accuracy
    metrics['precision'] = precision_score(y_true, y_pred, zero_division=0)
    metrics['recall'] = recall_score(y_true, y_pred, zero_division=0)
    metrics['f1_score'] = f1_score(y_true, y_pred, zero_division=0)
    
    # Probability-based metrics
    if len(np.unique(y_true)) == 2:  # Only for binary classification
        metrics['roc_auc'] = roc_auc_score(y_true, y_probs)
        metrics['pr_auc'] = average_precision_score(y_true, y_probs)
    
    # Confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    metrics['confusion_matrix'] = cm
    
    if cm.shape == (2, 2):
        tn, fp, fn, tp = cm.ravel()
        metrics['true_negatives'] = tn
        metrics['false_positives'] = fp
        metrics['false_negatives'] = fn
        metrics['true_positives'] = tp
        metrics['specificity'] = tn / (tn + fp) if (tn + fp) > 0 else 0
        metrics['sensitivity'] = tp / (tp + fn) if (tp + fn) > 0 else 0
    
    return metrics


def plot_results(df: pd.DataFrame, metrics: dict, output_dir: str = "vitis_cre/src/results/model_eval/self_predict"):
    """
    Create visualizations of model performance.
    
    Args:
        df: DataFrame with predictions
        metrics: Dictionary of calculated metrics
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Set up the plotting style with custom colors
    plt.style.use('default')
    sns.set_palette(CUSTOM_COLORS)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Add main title with validation methodology explanation
    fig.suptitle('DeepCRE Self-Prediction Performance Assessment\nEach model tested on its own validation chromosome (non-homologous gene set)', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 1. Confusion Matrix
    ax1 = axes[0, 0]
    cm = metrics['confusion_matrix']
    # Create custom colormap using our colors
    cmap = sns.blend_palette([CUSTOM_COLORS[0], CUSTOM_COLORS[1]], as_cmap=True)
    sns.heatmap(cm, annot=True, fmt='d', cmap=cmap, ax=ax1, 
                xticklabels=['Predicted 0', 'Predicted 1'],
                yticklabels=['True 0', 'True 1'])
    ax1.set_title('Confusion Matrix')
    ax1.set_ylabel('True Labels')
    ax1.set_xlabel('Predicted Labels')
    
    # 2. Probability Distribution by Class
    ax2 = axes[0, 1]
    for i, class_label in enumerate([0, 1]):
        class_probs = df[df['true_targets'] == class_label]['pred_probs']
        ax2.hist(class_probs, bins=30, alpha=0.7, 
                label=f'True Class {class_label}', density=True, 
                color=CUSTOM_COLORS[i])
    ax2.axvline(x=0.5, color=CUSTOM_COLORS[3], linestyle='--', linewidth=2, label='Threshold (0.5)')
    ax2.set_xlabel('Predicted Probability')
    ax2.set_ylabel('Density')
    ax2.set_title('Probability Distribution by True Class')
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
    
    # 4. Metrics Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Create metrics text with both decimal and percentage accuracy
    metrics_text = f"""
    PERFORMANCE METRICS
    
    Overall Accuracy: {metrics['accuracy']:.3f} ({metrics['accuracy_percent']:.1f}%)
    Precision: {metrics['precision']:.3f}
    Recall (Sensitivity): {metrics['recall']:.3f}
    Specificity: {metrics.get('specificity', 0):.3f}
    F1-Score: {metrics['f1_score']:.3f}
    
    ROC-AUC: {metrics.get('roc_auc', 0):.3f}
    PR-AUC: {metrics.get('pr_auc', 0):.3f}
    
    True Positives: {metrics.get('true_positives', 0):,}
    True Negatives: {metrics.get('true_negatives', 0):,}
    False Positives: {metrics.get('false_positives', 0):,}
    False Negatives: {metrics.get('false_negatives', 0):,}
    """
    
    # White background with subtle border
    ax4.text(0.1, 0.9, metrics_text, transform=ax4.transAxes, fontsize=12,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white", 
                      edgecolor=CUSTOM_COLORS[0], linewidth=1.5, alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/model_performance_summary.png", dpi=300, bbox_inches='tight')
    plt.show()


def threshold_analysis(df: pd.DataFrame, output_dir: str = "vitis_cre/src/results/model_eval/self_predict"):
    """
    Analyze performance across different probability thresholds.
    
    Args:
        df: DataFrame with predictions
        output_dir: Directory to save plots
    """
    thresholds = np.arange(0.1, 0.95, 0.05)
    threshold_metrics = []
    
    for threshold in thresholds:
        y_pred_thresh = (df['pred_probs'] >= threshold).astype(int)
        
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
    
    # Plot threshold analysis with custom colors
    plt.figure(figsize=(12, 8))
    
    # Add figure title explaining validation methodology
    plt.suptitle('DeepCRE Self-Prediction Threshold Analysis\nPerformance on validation chromosomes (non-homologous genes)', 
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
    plt.title('Model Performance vs Probability Threshold')
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
    plt.savefig(f"{output_dir}/threshold_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Find optimal threshold (max F1-score)
    optimal_idx = threshold_df['f1_score'].idxmax()
    optimal_threshold = threshold_df.loc[optimal_idx, 'threshold']
    optimal_f1 = threshold_df.loc[optimal_idx, 'f1_score']
    
    print(f"\nOptimal threshold for F1-score: {optimal_threshold:.2f} (F1 = {optimal_f1:.3f})")
    
    return threshold_df


def main():
    parser = argparse.ArgumentParser(description='Assess DeepCRE model accuracy')
    parser.add_argument('--predictions', '-p', required=True, 
                       help='Path to predictions CSV file')
    parser.add_argument('--threshold', '-t', type=float, default=0.5,
                       help='Probability threshold for binary classification (default: 0.5)')
    parser.add_argument('--output', '-o', default='vitis_cre/src/results/model_eval/self_predict',
                       help='Output directory for plots and results (default: vitis_cre/src/results/model_eval/self_predict)')
    
    args = parser.parse_args()
    
    # Load predictions
    print("Loading predictions...")
    df = load_predictions(args.predictions, args.threshold)
    
    # Calculate metrics
    print("\nCalculating metrics...")
    metrics = calculate_metrics(df['true_targets'], df['pred_binary'], df['pred_probs'])
    
    # Print summary
    print("\n" + "="*50)
    print("MODEL PERFORMANCE SUMMARY")
    print("="*50)
    print(f"Overall Accuracy: {metrics['accuracy']:.3f} ({metrics['accuracy_percent']:.1f}%)")
    print(f"Precision: {metrics['precision']:.3f}")
    print(f"Recall: {metrics['recall']:.3f}")
    print(f"F1-Score: {metrics['f1_score']:.3f}")
    print(f"ROC-AUC: {metrics.get('roc_auc', 'N/A'):.3f}")
    print(f"PR-AUC: {metrics.get('pr_auc', 'N/A'):.3f}")
    print("="*50)
    
    # Create visualizations
    print("\nGenerating plots...")
    plot_results(df, metrics, args.output)
    
    # Threshold analysis
    print("\nPerforming threshold analysis...")
    threshold_analysis(df, args.output)
    
    # Save detailed classification report
    report = classification_report(df['true_targets'], df['pred_binary'], 
                                 target_names=['Low Expression', 'High Expression'])
    
    with open(f"{args.output}/classification_report.txt", 'w') as f:
        f.write("DeepCRE Model Self-Prediction Classification Report\n")
        f.write("Each model tested on its own validation chromosome\n")
        f.write("(Non-homologous gene set)\n")
        f.write("="*50 + "\n\n")
        f.write(report)
    
    print(f"\nResults saved to {args.output}/")
    print("Files created:")
    print("- model_performance_summary.png")
    print("- threshold_analysis.png") 
    print("- classification_report.txt")


if __name__ == "__main__":
    main()