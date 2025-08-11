'''
This is an optional script (written after DeepCRE reimplementation) to extract motifs from modisco results.
It generates visualizations of motifs and saves them in a structured format.
Moca_blue pipeline can run without it
'''

import os
import h5py 
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from modisco.visualization.viz_sequence import plot_weights 
from modisco.visualization.viz_sequence import plot_weights_given_ax


script_dir = os.path.dirname(os.path.abspath(__file__))
modisco_dir = os.path.join(script_dir, 'results/modisco')
os.chdir(modisco_dir)

if not os.path.exists('figures/motifs/'):
    os.mkdir('figures/motifs/')


# ----------------- additional function for safe plotting and logging ------------------ #

def safe_plot_weights(matrix, save_path, pattern_name, strand, skipped_patterns):
    """
    Only plot if matrix is 2D and non-empty. Log any skipped patterns.
    """
    try:
        if matrix.ndim != 2 or matrix.size == 0:
            raise ValueError("Matrix not 2D or empty")

        fig = plt.figure(figsize=(20, 2))
        ax = fig.add_subplot(111)
        plot_weights_given_ax(ax=ax, array=matrix)
        plt.savefig(save_path, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        skipped_patterns.append({
            'pattern': pattern_name,
            'strand': strand,
            'shape': matrix.shape,
            'error': str(e),
            'save_path': save_path
        })


# ------------- Extracting motifs from modisco results for visualisation ---------------- #
def get_predictive_pwms(mod_file, specie):
    print(specie)
    cwm = []
    motif_id = []
    strand = []
    metacluster_id = []
    n_motif_seqlets = []

    skipped_patterns = []

    f = h5py.File(mod_file, "r")
    for metacluster_idx, metacluster_key in enumerate(f["metacluster_idx_to_submetacluster_results"]):
        metacluster = f["metacluster_idx_to_submetacluster_results"][metacluster_key]
        patterns = metacluster['seqlets_to_patterns_result']['patterns']
        print(metacluster_idx, metacluster_key, len(patterns['all_pattern_names']))

        for pattern_idx, pattern_name in enumerate(patterns['all_pattern_names']):
            pattern = patterns[pattern_name.decode()]
            pattern_seqlets = pattern["seqlets_and_alnmts"]["seqlets"]
            motif_name = pattern_name.decode()

            # Forward strand
            nfcwm = np.absolute(pattern["task0_contrib_scores"]["fwd"][:].astype(np.float32))
            nfcwm = len(pattern_seqlets) * (nfcwm / np.max(nfcwm.flat))
            cwm.append(nfcwm)
            motif_id.append(motif_name)
            metacluster_id.append(metacluster_key)
            n_motif_seqlets.append(len(pattern_seqlets))
            strand.append('fwd')
            save_fwd = f'figures/motifs/{specie}/{metacluster_key}_{motif_name}_fwd.png'
            safe_plot_weights(pattern["task0_contrib_scores"]["fwd"][:], save_fwd, motif_name, 'fwd', skipped_patterns)

            # Reverse strand
            nrcwm = np.absolute(pattern["task0_contrib_scores"]["rev"][:].astype(np.float32))
            nrcwm = len(pattern_seqlets) * (nrcwm / np.max(nrcwm.flat))
            cwm.append(nrcwm)
            motif_id.append(motif_name)
            metacluster_id.append(metacluster_key)
            n_motif_seqlets.append(len(pattern_seqlets))
            strand.append('rev')
            save_rev = f'figures/motifs/{specie}/{metacluster_key}_{motif_name}_rev.png'
            safe_plot_weights(pattern["task0_contrib_scores"]["rev"][:], save_rev, motif_name, 'rev', skipped_patterns)

    # Save cwms in H5 file
    cwm = np.array(cwm)
    motif_save = f'figures/motifs/{specie}/motifs.h5'
    if os.path.exists(motif_save):
        os.remove(motif_save)
    with h5py.File(motif_save, 'w') as h:
        h.create_dataset('CWMs', data=cwm)

    meta_info = pd.DataFrame([motif_id, strand, metacluster_id, n_motif_seqlets]).T
    meta_info.columns = ['motifID', 'strand', 'metacluster', 'n_seqlets']
    meta_info.to_csv(f'figures/motifs/{specie}/meta_info.csv', sep='\t', index=False)

# -------------------- Run on a single file --------------------- #
modisco_feats = ['vitis_02_modisco.hdf5']
species = ['VT']

for feats, sp in zip(modisco_feats, species):
    if not os.path.exists(f'figures/motifs/{sp}'):
        os.mkdir(f'figures/motifs/{sp}')
    feats_path = f'{feats}'
    get_predictive_pwms(feats_path, sp)
