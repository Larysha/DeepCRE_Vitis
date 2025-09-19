# script to extract CWMs from a MODISCO HDF5 file and compare them with JASPAR motifs using Tomtom
# requires h5py, numpy, and subprocess for running Docker commands (memesuite/memesuite)
# run after modisco output has been generated - before the moca_blue pipeline

import os
import h5py
import numpy as np
import subprocess
import glob

def find_modisco_file(results_dir="results/modisco"):
    """Find the most recent MoDISco HDF5 file"""
    pattern = os.path.join(results_dir, "*_modisco.hdf5")
    modisco_files = glob.glob(pattern)
    
    if not modisco_files:
        raise FileNotFoundError(f"No MoDISco HDF5 files found in {results_dir}")
    
    # Return the most recent file if multiple exist
    modisco_file = max(modisco_files, key=os.path.getmtime)
    print(f"[info] Using MoDISco file: {modisco_file}")
    return modisco_file

def find_jaspar_db(jaspar_dir="jaspar"):
    """Find JASPAR combined MEME database"""
    possible_names = [
        "jaspar_plants_2024_combined.meme",
        "jaspar_core_plants_2024_combined.meme",
        "jaspar_combined.meme"
    ]
    
    for name in possible_names:
        jaspar_path = os.path.join(jaspar_dir, name)
        if os.path.exists(jaspar_path):
            print(f"[info] Using JASPAR database: {jaspar_path}")
            return jaspar_path
    
    raise FileNotFoundError(f"No JASPAR database found in {jaspar_dir}. Expected one of: {possible_names}")

def extract_all_cwms_from_modisco(h5_path, min_len=4):
    """Extract all CWM matrices from MoDISco HDF5 file"""
    if not os.path.exists(h5_path):
        raise FileNotFoundError(f"MoDISco file not found: {h5_path}")
    
    cwms, names = [], []
    with h5py.File(h5_path, "r") as f:
        mc_keys = list(f["metacluster_idx_to_submetacluster_results"].keys())
        print(f"[info] Processing {len(mc_keys)} metaclusters: {mc_keys}")
        
        for mc_key in mc_keys:
            patterns = f["metacluster_idx_to_submetacluster_results"][mc_key]["seqlets_to_patterns_result"]["patterns"]

            for pattern_key in sorted(patterns.keys()):
                if not pattern_key.startswith("pattern_"):
                    continue

                try:
                    cwm = patterns[pattern_key]["sequence"]["modisco_cwm"][:]
                    source = f"{mc_key}/modisco_cwm"
                except KeyError:
                    try:
                        cwm = patterns[pattern_key]["task0_contrib_scores"]["fwd"][:]
                        source = f"{mc_key}/task0_contrib_scores/fwd"
                    except KeyError:
                        print(f"[warn] No usable matrix for {mc_key}/{pattern_key}, skipping.")
                        continue

                cwm = np.array(cwm)

                if cwm.ndim > 2:
                    if 4 in cwm.shape:
                        axis_to_keep = np.where(np.array(cwm.shape) == 4)[0][0]
                        cwm = np.moveaxis(cwm, axis_to_keep, -1)
                        cwm = cwm.reshape(-1, 4)
                    else:
                        print(f"[warn] {mc_key}/{pattern_key}: unexpected shape {cwm.shape}, skipping.")
                        continue
                elif cwm.shape[0] == 4 and cwm.shape[1] != 4:
                    cwm = cwm.T

                if cwm.shape[1] != 4 or cwm.shape[0] < min_len:
                    print(f"[warn] Skipping {mc_key}/{pattern_key} due to invalid shape {cwm.shape}")
                    continue

                motif_name = f"{mc_key}_{pattern_key}"
                print(f"[info] {motif_name}: shape {cwm.shape} from {source}")
                cwms.append(cwm)
                names.append(motif_name)

    print(f"[info] Extracted {len(cwms)} valid CWMs across all metaclusters")
    return cwms, names


def write_cwms_to_meme_file(cwms, names, output_file):
    """Write CWM matrices to MEME format file"""
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, "w") as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies:\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        for name, cwm in zip(names, cwms):
            if cwm.shape[0] == 4 and cwm.shape[1] != 4:
                cwm = cwm.T

            motif_len = cwm.shape[0]
            f.write(f"MOTIF {name}\n")
            f.write(f"letter-probability matrix: alength= 4 w= {motif_len} nsites= 20 E= 0\n")

            for row in cwm:
                row = np.clip(row, 0, None)
                row_sum = np.sum(row)
                if row_sum == 0:
                    row = np.array([0.25, 0.25, 0.25, 0.25])
                else:
                    row = row / row_sum
                f.write(" ".join(f"{x:.6f}" for x in row) + "\n")
            f.write("\n")

    print(f"[info] Wrote MEME motif file to: {output_file}")


def run_tomtom_docker(query_meme, target_meme, output_dir):
    """Run TOMTOM comparison via Docker"""
    # Check if Docker is available
    try:
        subprocess.run(["docker", "--version"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("Docker is not available. Install Docker or run TOMTOM manually.")
    
    # Get absolute paths
    abs_query = os.path.abspath(query_meme)
    abs_target = os.path.abspath(target_meme)
    abs_out = os.path.abspath(output_dir)
    
    # Find a common parent directory to mount
    # Use the current working directory as the mount point
    mount_dir = os.getcwd()
    
    # Calculate relative paths from mount directory
    rel_query = os.path.relpath(abs_query, mount_dir)
    rel_target = os.path.relpath(abs_target, mount_dir)
    rel_output = os.path.relpath(abs_out, mount_dir)

    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{mount_dir}:/data",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "memesuite/memesuite",
        "tomtom",
        "-oc", f"/data/{rel_output}",
        f"/data/{rel_query}",
        f"/data/{rel_target}"
    ]

    print(f"[info] Running TOMTOM via Docker...")
    print(f"[info] Command: {' '.join(docker_cmd)}")
    
    try:
        subprocess.run(docker_cmd, check=True)
        print(f"[info] TOMTOM results written to: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"[error] TOMTOM failed with return code {e.returncode}")
        raise


def main():
    """Main execution function"""
    print("[info] Starting MoDISco to JASPAR comparison pipeline...")
    
    # Define paths relative to project root (~/phd_2025/deepcre_vitis/vitis_cre/src)
    try:
        modisco_h5 = find_modisco_file("results/modisco")
        jaspar_db = find_jaspar_db("jaspar")
        
        output_dir = "results/tomtom"
        extracted_meme = os.path.join(output_dir, "extracted_cwms.meme")
        
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"[info] Output directory: {output_dir}")
        
        # Execute pipeline
        cwms, names = extract_all_cwms_from_modisco(modisco_h5)
        
        if len(cwms) == 0:
            print("[warn] No motifs extracted. Check your MoDISco file.")
            return
        
        write_cwms_to_meme_file(cwms, names, extracted_meme)
        run_tomtom_docker(extracted_meme, jaspar_db, output_dir)
        
        print("\n" + "="*50)
        print("COMPARISON COMPLETE")
        print("="*50)
        print(f"Results saved to: {output_dir}")
        print("Next steps:")
        print("1. Check tomtom.html for interactive results")
        print("2. Run annotate_tomtom.py to add TF metadata")
        print("3. Proceed with moca_blue pipeline")
        
    except Exception as e:
        print(f"[error] Pipeline failed: {e}")
        print("\nTroubleshooting:")
        print("- Ensure MoDISco HDF5 file exists in results/modisco/")
        print("- Ensure JASPAR database exists in jaspar/")
        print("- Check Docker is installed and running")
        raise


if __name__ == "__main__":
    main()