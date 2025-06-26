import os
import h5py
import numpy as np
import subprocess

def extract_all_cwms_from_modisco(h5_path, min_len=4):
    cwms, names = [], []
    with h5py.File(h5_path, "r") as f:
        mc_keys = list(f["metacluster_idx_to_submetacluster_results"].keys())
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
    abs_query = os.path.abspath(query_meme)
    abs_target = os.path.abspath(target_meme)
    abs_out = os.path.abspath(output_dir)
    modisco_dir = os.path.dirname(abs_target)

    rel_query = os.path.relpath(abs_query, modisco_dir)
    rel_target = os.path.relpath(abs_target, modisco_dir)
    rel_output = os.path.relpath(abs_out, modisco_dir)

    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{modisco_dir}:/data",
        "-u", f"{os.getuid()}:{os.getgid()}",
        "memesuite/memesuite",
        "tomtom",
        "-oc", f"/data/{rel_output}",
        f"/data/{rel_query}",
        f"/data/{rel_target}"
    ]

    print(f"[info] Running Tomtom via Docker:\n{' '.join(docker_cmd)}")
    subprocess.run(docker_cmd, check=True)
    print(f"[info] Tomtom results written to: {output_dir}")


if __name__ == "__main__":
    modisco_h5 = "../results/modisco/vitis_02_modisco.hdf5"
    output_dir = "../results/modisco/tomtom_results"
    extracted_meme = os.path.join(output_dir, "extracted_cwms.meme")
    jaspar_db = "../results/modisco/jaspar_plants_2024_combined.meme"

    os.makedirs(output_dir, exist_ok=True)

    cwms, names = extract_all_cwms_from_modisco(modisco_h5)
    write_cwms_to_meme_file(cwms, names, extracted_meme)
    run_tomtom_docker(extracted_meme, jaspar_db, output_dir)
