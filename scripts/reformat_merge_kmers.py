import argparse
from pathlib import Path

import polars as pl  # type: ignore
from tqdm import tqdm  # type: ignore
import pandas as pd  # type: ignore

def reverse_complement(dna_sequence: str) -> str:
    """
    Compute the reverse complement of a DNA sequence in upper case.

    Parameters:
    dna_sequence: DNA sequence in upper case (only ATCG)

    Returns:
    Reverse complement of dna_sequence.
    """
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
    }
    complement_sequence = "".join(complement[base] for base in dna_sequence)
    return complement_sequence[::-1]


def derive_sample_name(file_path: str) -> str:
    """
    Keep backward-compatible sample naming from *.kmer.tsv style paths.
    """
    path = Path(file_path)
    if path.name.endswith(".kmer.tsv"):
        return path.name[: -len(".kmer.tsv")]
    return path.stem


def load_kmer_table(file_path: str, sample_name: str) -> pl.DataFrame:
    with open(file_path, "r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()

    if first_line == f"KMER\t{sample_name}":
        kmer_df = pl.read_csv(file_path, separator="\t")
    else:
        kmer_df = pl.read_csv(
            file_path,
            separator=" ",
            has_header=False,
            new_columns=["KMER", sample_name],
        )

    return kmer_df.with_columns(pl.col(sample_name).cast(pl.Float64, strict=False))


def main() -> None:
    
    ### manage parameters
    parser = argparse.ArgumentParser(
        description=(
            "Read k-mer tables from multiple samples, reformat each table and merges table across samples"
        )
    )
    parser.add_argument(
        "--kmer_file_list",
        required=True,
        help="Path to a text file containing one kmer file path per line",
    )
    parser.add_argument(
        "--kmer_info_file",
        required=True,
        help="Path to a file containing kmer cluster information (e.g., data/all_tagginer_kmers.tsv.gz)",
    )
    args = parser.parse_args()

    
    ### manage inputs
    kmer_files: list[dict[str, str]] = []
    with open(args.kmer_file_list, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            kmer_files.append({"NAME": derive_sample_name(line), "FILE": line})


    ### read and reformat kmer tables
    reformat_kmer_tables = []
    for row in tqdm(kmer_files, total=len(kmer_files), desc="Loading and reformatting k-mer files"):
        ### load kmer table
        kmer_df = load_kmer_table(row["FILE"], row["NAME"])

        ### make table with reverse complement
        rc_kmer_df = kmer_df.with_columns(
            pl.col("KMER").map_elements(reverse_complement, return_dtype=pl.String)
        )

        ### merge original and reverse complement tables, keeping unique kmers only
        kmer_df_expanded = (
            pl.concat([kmer_df, rc_kmer_df], how="vertical")
            .unique(keep="first")
            .sort("KMER")
        )

        ### append sample kmers to list
        reformat_kmer_tables.append(kmer_df_expanded)
        
    ### merge kmer tables across samples
    if not reformat_kmer_tables:
        raise ValueError("No input k-mer files found in --kmer_file_list.")
    # check that all tables have the same kmers
    reference_kmers = reformat_kmer_tables[0].get_column("KMER")
    for idx, table in enumerate(reformat_kmer_tables[1:], start=1):
        if not table.get_column("KMER").equals(reference_kmers):
            raise ValueError(
                f"KMER order/content mismatch for table index {idx}. "
                "All tables must have identical sorted KMER values."
            )
    # merge into a single table
    merged_kmer_df = pl.concat(
        [reformat_kmer_tables[0].select("KMER")]
        + [table.select(table.columns[1]) for table in reformat_kmer_tables],
        how="horizontal",
    )

    ### load data/all_tagging_kmers.tsv.gz
    info_df = pd.read_csv(args.kmer_info_file, sep="\t")[["KMER", "CLUSTER"]].set_index("KMER")
    # merge counts with information
    merged = info_df.join(merged_kmer_df.to_pandas().set_index("KMER"), how="left")
    
    merged.to_csv("final_kmer_merged.tsv", sep="\t", index=True)    

if __name__ == "__main__":
    main()
