# script to fetch metadata for JASPAR matrices and annotate a TOMTOM results file.
# load `tomtom.tsv` - filters significat matches (only E-value ≤ 0.05 AND q-value ≤ 0.05)
# calls JASPAR API to fetch TF names, families, and species
# combines TOMTOM scores with TF annotations
# exports annotated results to `tomtom_annotated.tsv`

import requests
import pandas as pd
import time
import os

JASPAR_API = "https://jaspar.genereg.net/api/v1/matrix/{}?format=json"

def fetch_jaspar_metadata(matrix_id):
    """Fetch TF name/species via JASPAR API; handles both dict and list responses."""
    url = JASPAR_API.format(matrix_id)
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()

        # DEBUG: print raw data returned by the API
        print(f"[debug] Raw response for {matrix_id}: {data}")

        # If list, take the first entry
        if isinstance(data, list):
            data = data[0]

        return {
            "Target_ID": matrix_id,
            "TF_name": data.get("name"),
            "TF_family": data.get("family", {}).get("name") if isinstance(data.get("family"), dict) else None,
            "Species": data.get("species", {}).get("name") if isinstance(data.get("species"), dict) else None
        }

    except Exception as e:
        print(f"[warn] Failed to fetch metadata for {matrix_id}: {e}")
        return {
            "Target_ID": matrix_id,
            "TF_name": None,
            "TF_family": None,
            "Species": None
        }


def main():
    # Define paths relative to project root
    tomtom_input = "results/tomtom/tomtom.tsv"
    tomtom_output = "results/tomtom/tomtom_annotated.tsv"
    
    # Check if input file exists
    if not os.path.exists(tomtom_input):
        print(f"[error] Input file not found: {tomtom_input}")
        print("Make sure you've run TOMTOM comparison first.")
        return
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(tomtom_output), exist_ok=True)
    
    print(f"[info] Reading TOMTOM results from: {tomtom_input}")
    df = pd.read_csv(tomtom_input, sep="\t", comment="#")
    
    # Filter for high-confidence matches
    df_f = df[(df["E-value"] <= 0.05) & (df["q-value"] <= 0.05)].copy()
    print(f"[info] {len(df_f)} high-confidence matches retained (E≤0.05 & q≤0.05) out of {len(df)} total matches.")
    
    if len(df_f) == 0:
        print("[warn] No significant matches found. Consider relaxing thresholds.")
        return
    
    unique_targets = df_f["Target_ID"].unique()
    print(f"[info] Fetching metadata for {len(unique_targets)} unique JASPAR matrices...")
    
    meta_list = []
    for i, mid in enumerate(unique_targets, 1):
        print(f"[info] Fetching {i}/{len(unique_targets)}: {mid}")
        meta_list.append(fetch_jaspar_metadata(mid))
        # Be nice to the JASPAR API
        time.sleep(0.5)
    
    meta_df = pd.DataFrame(meta_list)
    
    df_out = df_f.merge(meta_df, on="Target_ID", how="left")
    df_out.to_csv(tomtom_output, sep="\t", index=False)
    print(f"[info] Annotated results written to: {tomtom_output}")
    
    # Print summary of annotations
    annotated_count = df_out["TF_name"].notna().sum()
    print(f"[info] Successfully annotated {annotated_count}/{len(df_out)} matches with TF information.")


if __name__ == "__main__":
    main()