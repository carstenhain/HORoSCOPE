import argparse
import re
import pandas as pd # type: ignore
import pickle
import numpy as np # type: ignore
from tqdm import tqdm # type: ignore
import warnings
warnings.filterwarnings(
    "ignore",
    message="X does not have valid feature names, but PCA was fitted with feature names"
)
warnings.filterwarnings(
    "ignore",
    message=r"Trying to unpickle estimator .*",
    module=r"sklearn\.base"
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kmer_table', required=True, help="Table with kmer counts for each sample, output of reformat_merge_kmers.py")
    parser.add_argument('--samplesheet', required=True, help="Samplesheet with sample names, paths and normalization strategy, same as for the complete workflow")
    parser.add_argument('--model_directory', required=True, help="Directory containing the trained models for length inference")
    args = parser.parse_args()

    ### process only selected chromosomes
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr16", "chr17", "chr18", "chr19", "chr20"]

    ### prepare list to store results
    results = []

    ### load kmer counts and kmer annotations
    kmer_df = pd.read_csv(args.kmer_table, sep="\t").set_index("KMER")
    
    ### load samplesheet
    samplesheet = pd.read_csv(args.samplesheet, sep=";")
    
    
    
    
    
    ####### NORMALIZATION
    
    
    ### calculate normalization
    sample_names = samplesheet["NAME"].tolist()
    norm_dfs = []
    
    ### factor for different normalization values between raw and dbg samples, empirically determined
    raw_to_dbg_factor = 1.0475077160727808
    
    ### calculate global normalization value
    global_norm_df = kmer_df[kmer_df["CLUSTER"].str.contains("_NORM_")][sample_names].agg(["mean", "std"]).T
    global_norm_df["NORM_STRATEGY"] = "GLOBAL"
    norm_dfs.append(global_norm_df.copy())
    
    ### calculate chromosome-specific normalization values
    for chrom in chromosomes:
        ### all normalization kmers of this chromosome
        chrom_norm_df = kmer_df[kmer_df["CLUSTER"].str.contains(f"{chrom}_NORM_")][sample_names].agg(["mean", "std"]).T
        chrom_norm_df["NORM_STRATEGY"] = chrom
        norm_dfs.append(chrom_norm_df.copy())
        
        ### only p or q normalization kmers of this chromosome
        for arm in ["P", "Q"]:           
            chrom_norm_df = kmer_df[kmer_df["CLUSTER"] == f"{chrom}_NORM_{arm}"][sample_names].agg(["mean", "std"]).T
            chrom_norm_df["NORM_STRATEGY"] = f"{chrom}_{arm}"
            norm_dfs.append(chrom_norm_df.copy())
    
    ### merge normalization values into a single dataframe
    norm_df = pd.concat(norm_dfs)
    
    ### adjust normalization values for raw samples
    col_normed_mean = []
    for idx, row in norm_df.iterrows():
        if samplesheet[samplesheet["NAME"] == idx]["FILE_TYPE"].values[0] in ["fastq", "fasta", "fastq.gz", "fasta.gz", "fastq.bgz", "fasta.bgz", "bam", "cram"]:
            col_normed_mean.append(raw_to_dbg_factor * row["mean"])
        else:
            col_normed_mean.append(row["mean"])
    norm_df["mean"] = col_normed_mean
    
    ### apply different normalization strategies to kmer counts
    # do so for each chromosomes separatly, concat afterwards
    normed_kmer_dfs = []
    for chrom in tqdm(chromosomes, total=len(chromosomes), desc="Normalizing k-mer counts by chromosome"):
        
        ### copy kmers of this chromosome
        chrom_kmers = kmer_df[kmer_df["CLUSTER"].str.startswith(f"{chrom}_")].copy()
        
        ### build normalization vector
        normalization_vector = []
        for sample_column in chrom_kmers.columns:
            if sample_column == "CLUSTER":
                continue
            
            norm_strategy = samplesheet[samplesheet["NAME"] == sample_column]["NORMALIZATION"].values[0]
            if norm_strategy == "global":
                normalization_value = norm_df[norm_df["NORM_STRATEGY"] == "GLOBAL"].at[sample_column, "mean"]
            elif norm_strategy == "chromosome":
                normalization_value = norm_df[norm_df["NORM_STRATEGY"] == chrom].at[sample_column, "mean"]
            elif norm_strategy == "p_arm":
                normalization_value = norm_df[norm_df["NORM_STRATEGY"] == f"{chrom}_P"].at[sample_column, "mean"]
            elif norm_strategy == "q_arm":
                normalization_value = norm_df[norm_df["NORM_STRATEGY"] == f"{chrom}_Q"].at[sample_column, "mean"]
            else:
                ### test for external normalization strategy
                try:
                    normalization_value = float(norm_strategy)
                except ValueError:
                    raise ValueError(f"Unknown or erroneous normalization strategy '{norm_strategy}' for sample '{sample_column}', must be one of 'global', 'chromosome', 'p_arm', 'q_arm' or a numeric value.")
            
            normalization_vector.append({"SAMPLE": sample_column, "VALUE": normalization_value})
        
        ### convert normalization vector into dataframe and apply to kmer counts
        chrom_norm_df = pd.DataFrame(normalization_vector).set_index("SAMPLE")
        chrom_kmers[chrom_norm_df.index] = chrom_kmers[chrom_norm_df.index].div(chrom_norm_df["VALUE"], axis=1)
        
        ### append normed kmers of this chromosome to list
        normed_kmer_dfs.append(chrom_kmers)
    
    ### concat normed kmer dataframes of all chromosomes
    normed_kmer_df = pd.concat(normed_kmer_dfs)

    ### build dataframe with normalization metrics for QC
    ### precompute kmers sets
    norm_kmer_sets = {
        "GLOBAL": normed_kmer_df[normed_kmer_df["CLUSTER"].str.contains("_NORM_")].index.tolist()
    }
    for chrom in chromosomes:
        norm_kmer_sets[chrom] = normed_kmer_df[normed_kmer_df["CLUSTER"].str.contains(f"{chrom}_NORM_")].index.tolist()
        norm_kmer_sets[f"{chrom}_P"] = normed_kmer_df[normed_kmer_df["CLUSTER"].str.contains(f"{chrom}_NORM_P")].index.tolist()
        norm_kmer_sets[f"{chrom}_Q"] = normed_kmer_df[normed_kmer_df["CLUSTER"].str.contains(f"{chrom}_NORM_Q")].index.tolist()
        
    
    ### precompute number of found kmers for each strategy and sample
    found_kmer_cutoff = 0
    strategy_counts = {
        strategy: (normed_kmer_df.loc[kmers, sample_names] > found_kmer_cutoff).sum(axis=0)
        for strategy, kmers in norm_kmer_sets.items()
    }

    ### add normalization metrics to norm_df
    col_num_norm_kmers_found = []
    col_frac_norm_kmers_found = []
    
    for idx, row in tqdm(norm_df.iterrows(), desc="Calculating normalization metrics", total=norm_df.shape[0]):

        ### number of found kmers are those with a count > 0
        strategy = row["NORM_STRATEGY"]
        num_found = strategy_counts[strategy].at[idx]
        
        col_num_norm_kmers_found.append(num_found)
        col_frac_norm_kmers_found.append(num_found / len(norm_kmer_sets[strategy]))
        
    norm_df["NUM_KMERS_FOUND"] = col_num_norm_kmers_found
    norm_df["FRAC_KMERS_FOUND"] = col_frac_norm_kmers_found
    
    ### write normalization metrics to file for QC
    norm_df[["NORM_STRATEGY", "mean", "std", "NUM_KMERS_FOUND", "FRAC_KMERS_FOUND"]].to_csv("normalization_metrics.tsv", sep="\t", index=True)
    
    
    
    
    
    ####### GENOTYPING
    
    # load prediction models
    kmer_sets = pickle.load(open(f"{args.model_directory}/length_kmers.pkl", "rb"))
    pca_models = pickle.load(open(f"{args.model_directory}/pca_models.pkl", "rb"))
    linreg_models = pickle.load(open(f"{args.model_directory}/linreg_models.pkl", "rb"))
    
    ### query all chromosomes and determine cluster, save the fraction of found kmers in a dataframe    
    for chrom in tqdm(chromosomes, desc="Genotyping by chromosome", total=len(chromosomes)):
        
        ### gather all clusters for this chromosomes
        pattern = re.compile(rf"^{re.escape(chrom)}_\d{{1,5}}$")
        cluster = [
            x
            for x in normed_kmer_df["CLUSTER"].dropna().astype(str).unique().tolist() # type: ignore
            if pattern.fullmatch(x)
        ]
        
        ### get fraction of found kmers (kmer count > 2) for each cluster
        for c in cluster:

            ### subset kmers to this cluster
            c_subset = normed_kmer_df[normed_kmer_df["CLUSTER"] == c][sample_names]
            
            ### get number of kmers and number of positive kmers
            c_kmers = c_subset.shape[0]

            ### skip clusters with very few tagging kmers
            if c_kmers < 10:
                c_positive = 0
            else:   
                c_positive = (c_subset > 0.2).sum()
            
            ### calculate fraction and add information about cluster and chromosome
            result_series = c_positive / c_kmers
            result_series = pd.Series(result_series)
            result_series["CLUSTER"] = c
            result_series["CHROM"] = chrom
            result_series["TYPE"] = "CLUSTER"

            ### append results
            results.append(result_series)
        
        
        ### predict length

        model_name = f"{chrom}_HOR"

        # skip if no model is present or not enough kmers are used for the model
        if not (model_name in kmer_sets):
            print(f"No model found for {chrom}, skipping length prediction.")
            continue
        if len(kmer_sets[model_name]) <= 100:
            print(f"Not enough kmers for {chrom} length prediction, skipping.")
            continue
        
        ### extract length kmers
        length_kmers = kmer_sets[model_name]

        ### prepare results
        result_series = pd.Series()
        result_series["CLUSTER"] = "LENGTH"
        result_series["CHROM"] = chrom
        result_series["TYPE"] = "LENGTH"

        ### predict length for each sample
        for sample in sample_names:
            sample_kmer_vector = normed_kmer_df[normed_kmer_df["CLUSTER"] == model_name].loc[length_kmers, sample].to_numpy(dtype=float).reshape(1, -1)
            total_length_reduced = pca_models[model_name].transform(sample_kmer_vector)
            total_length_prediction = linreg_models[model_name].predict(total_length_reduced)[0]
            result_series[sample] = total_length_prediction

        ### append results
        results.append(result_series)
        
    results_df = pd.DataFrame(results)

    #results_df.to_csv("centromere_genotyping_results.tsv", sep="\t", index=False)

    formatted_results = []

    for sample in sample_names:

        for chrom in chromosomes:

            cluster_passing = results_df[
                (results_df["TYPE"] == "CLUSTER") & 
                (results_df["CHROM"] == chrom) & 
                (results_df[sample] > 0.5)
            ]["CLUSTER"].tolist()

            
            if len(cluster_passing) == 0 or len(cluster_passing) > 2:
                h1 = "NA"
                h2 = "NA"
            if len(cluster_passing) == 1:
                h1 = cluster_passing[0]
                h2 = cluster_passing[0]
            if len(cluster_passing) == 2:
                h1 = sorted(cluster_passing)[0]
                h2 = sorted(cluster_passing)[1]

            formatted_results.append({
                "SAMPLE":sample,
                "CHROM":chrom,
                "CLUSTER_H1":h1,
                "CLUSTER_H2":h2,
                "HOR_LENGTH":results_df[(results_df["CHROM"] == chrom) & (results_df["CLUSTER"] == "LENGTH")][sample].values[0]
            })

    ### write to file
    pd.DataFrame(formatted_results).to_csv("centromere_genotyping.tsv", sep="\t", index=False)
    
if __name__ == '__main__':
    main()
