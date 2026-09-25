import pandas as pd
import numpy as np
from tqdm import tqdm
import os

def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    
    # Avoid division by zero if a vector is all zeros
    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0
    
    return dot_product / (norm_v1 * norm_v2)

def analyze_mutational_profiles(file_path, n_simulations=1000):
    print("Loading data...")
    # Load the dataset
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"Error: The file was not found at {file_path}")
        return None
    df_mutations = df.set_index("MutationType")

    # Store results for final DataFrame construction
    results = []
    
    # Use tqdm for a progress bar during the simulation
    print(f"Running {n_simulations} simulations for each of the {df_mutations.shape[1]} samples...")
    for sample_id in tqdm(df_mutations.columns, desc="Processing Samples"):
        original_profile = df_mutations[sample_id].values
        original_mutation_count = original_profile.sum()

        # Skip samples with no mutations
        if original_mutation_count == 0:
            print('Skipping sample with zero mutations:', sample_id)
            continue

        for i in range(n_simulations):
            # Add Poisson noise
            simulated_profile = np.random.poisson(original_profile)
            simulated_mutation_count = simulated_profile.sum()
            
            # Calculate cosine similarity
            cosim = cosine_similarity(original_profile, simulated_profile)
            
            # Append results
            results.append({
                "simulation_id": f"{sample_id}_{i}",
                "sample_id": sample_id,
                "original_mutation_count": original_mutation_count,
                "poisson_mutation_count": simulated_mutation_count,
                "cosine_similarity": cosim
            })

    print("Simulations complete.")
    # Create and return the final results DataFrame
    results_df = pd.DataFrame(results)
    return results_df


if __name__ == '__main__':
    file_path = 'combined_esignatures_samples.SBS96.all'
    output_filename = 'combined_esignatures_samples_1386_SBS96_results.csv'
    
    # Run the analysis
    results_df = analyze_mutational_profiles(file_path, n_simulations=1000)

    if results_df is not None:
        # Save Results 
        output_path = os.path.join(os.path.dirname(file_path), output_filename)
        print(f"Saving results to {output_path}...")
        print('Analyze and plot using R')
        results_df.to_csv(output_path, index=False)

        # Metric Calculations 
        print("Rough summary of results...")
        # 1. Calculate cosine similarity quantiles for the entire dataset
        quantiles = results_df['cosine_similarity'].quantile([0.25, 0.50, 0.75])
        print("\n Cosine Similarity Statistics ")
        print(f"\t25th Percentile: {quantiles[0.25]:.4f}")
        print(f"\t50th Percentile (Median): {quantiles[0.50]:.4f}")
        print(f"\t75th Percentile: {quantiles[0.75]:.4f}")

        # 2. Calculate mutation count thresholds
        print("\nCalculating mutation count thresholds...")
        
        # Create a boolean column for stability
        results_df['is_stable'] = results_df['cosine_similarity'] >= 0.95
        
        # Sort by simulated mutation count in descending order
        df_sorted = results_df.sort_values('poisson_mutation_count', ascending=False)
        
        # Calculate the cumulative number of stable simulations
        df_sorted['stable_cumsum'] = df_sorted['is_stable'].cumsum()
        
        # Calculate the cumulative total number of simulations
        df_sorted['total_cumsum'] = np.arange(1, len(df_sorted) + 1)
        
        # Calculate the proportion of stable simulations at or above each count
        df_sorted['stable_proportion'] = df_sorted['stable_cumsum'] / df_sorted['total_cumsum']
        
        # Find the threshold for 90%
        threshold_90_rows = df_sorted[df_sorted['stable_proportion'] >= 0.90]
        threshold_90 = threshold_90_rows['poisson_mutation_count'].min() if not threshold_90_rows.empty else "Not Found"
        
        # Find the threshold for 95%
        threshold_95_rows = df_sorted[df_sorted['stable_proportion'] >= 0.95]
        threshold_95 = threshold_95_rows['poisson_mutation_count'].min() if not threshold_95_rows.empty else "Not Found"

        print("\n Mutation Count Thresholds (for cos(sim) >= 0.95) ")
        print(f"\tThreshold for >= 90% of simulations to be stable: {threshold_90} mutations")
        print(f"\tThreshold for >= 95% of simulations to be stable: {threshold_95} mutations")
        print("\nRough summary complete. Analyze and plot it using R.")
    else:
        print("No results to save or analyze due to an error in loading results_df.")