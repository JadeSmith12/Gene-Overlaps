# Gene-Overlaps
Identifying gene overlaps with regions of interest ie. methylated/acetylated regions

import pandas as pd

# Read in genes of interest from CSV with format Chromosome, start, stop, protid
genes = pd.read_csv("Gene_Loc.csv")
​
# Read in methylated/acetylated regions from CSV with format Chromosome, start, stop
methylated_regions = pd.read_csv("Methylated_Regions.csv")

# Function to check for overlap between two regions
def check_overlap(chromosome1, start1, end1, chromosome2, start2, end2):
    return chromosome1 == chromosome2 and start1 <= end2 and end1 >= start2
    
# Initialize a counter for overlapping genes
overlap_count = 0

# Initialize a list to store the protids of overlapping genes
overlapping_protids = []

# Loop through each gene
for _, gene in genes.iterrows():
    # Extract chromosome, start, and end positions of the current gene
    chromosome_gene = gene["Chromosome"]
    start_gene = gene["start"]
    end_gene = gene["stop"]
    protid_gene = gene["protid"]
    
    # Check for overlap with each methylated region
    overlaps = methylated_regions.apply(lambda row: check_overlap(chromosome_gene, start_gene, end_gene, row["Chromosome"], row["start"], row["stop"]), axis=1)
    
    # If any overlap found, increment the counter and add the protid to the list
    if overlaps.any():
        overlap_count += 1
        overlapping_protids.append(protid_gene)
​
# Calculate the percentage of genes with overlap
total_genes = len(genes)
overlap_percentage = (overlap_count / total_genes) * 100

# Print the percentage of genes overlapping with methylated/acetylated regions
print("Percentage of genes overlapping with methylated regions:", overlap_percentage)


# Now generate a csv containing protein id and whether there is an overlap or not

# Initialize lists to store protids and corresponding overlap status
protids = []
overlap_status = []
​
# Loop through each gene
for _, gene in genes.iterrows():
    # Extract chromosome, start, and end positions of the current gene
    chromosome_gene = gene["Chromosome"]
    start_gene = gene["start"]
    end_gene = gene["stop"]
    protid_gene = gene["protid"]
    
    # Check for overlap with each methylated region
    overlaps = methylated_regions.apply(lambda row: check_overlap(chromosome_gene, start_gene, end_gene, row["Chromosome"], row["start"], row["stop"]), axis=1)
    
    # Check if any overlap found
    if overlaps.any():
        overlap = "Overlap"
    else:
        overlap = "No overlap"
    
    # Append protid and overlap status to the respective lists
    protids.append(protid_gene)
    overlap_status.append(overlap)
​
# Create a DataFrame to store the results
results_df = pd.DataFrame({"Protid": protids, "Overlap Status": overlap_status})
​
# Save the results to a CSV file
results_df.to_csv("overlap_results.csv", index=False)
