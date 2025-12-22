#!/usr/bin/env python3
# RBIF100 Week 7/8 - Final Project: Bioinformatics Data Analysis and Visualization
# Name: Helena Balbat

# Step-by-Step Tasks:
    # STEP 1 = retrieve biological data
    # STEP 2 - data processing & analysis
    # STEP 3 - data visualization
    # STEP 4 - summary report (see analysis_report.md)


# STEP 0 - import all necessary packages (throughout all steps)
import requests
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# STEP 1 = retrieve biological data
# Tasks:
    # Use a RESTful API to retrieve biological data
    # Select gene of interest or organism, but retrieve at least 100 data entries
        # NOTE: selected APOE (apolipoprotein E) as gene of interest
        # NOTE: retrieving genetic variants in Homo sapiens from NCBI
    # Save the data in an appropriate format (CSV)

# search parameters
gene_name = "APOE"
esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

params = {
    "db": "gene",  # search NCBI Gene database
    "term": f"{gene_name}[Gene Name] AND Homo sapiens[Organism]",  # human APOE
    "retmode": "json"
}

# send GET request to NCBI API
response = requests.get(esearch_url, params=params)
data = response.json()

# extract gene ID automatically
gene_id_list = data.get("esearchresult", {}).get("idlist", [])
if not gene_id_list:    # check if gene ID list is empty
    raise ValueError(f"No Gene ID found for {gene_name}")   # raise an error if no Gene ID is found for APOE
gene_id = gene_id_list[0]   # take the first gene ID (typical APOE gene)
print("\nStep 1: Retrieval\nNCBI Gene ID for APOE:", gene_id)    # print first gene ID to show successful data request

# use Entrez esearch to get dbSNP IDs linked to this gene (SNP database)
esearch_snp_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

params_snp = {
    "db": "snp",
    "term": f"APOE[GENE] AND Homo sapiens[ORGN]",  # use dynamically retrieved gene_id
    "retmax": 300,               # number of variants to retrieve (at least 100)
    "retmode": "json"
}

# send another GET request this time to dbSNP, store data
response_snp = requests.get(esearch_snp_url, params=params_snp)
data_snp = response_snp.json()

variant_ids = data_snp.get("esearchresult", {}).get("idlist", [])
print("Number of variants retrieved:", len(variant_ids))    # print how many variants were retrieved (should be <= max, which is 300 here)

# empty list to store processed variant records
variant_records = []

# loop over variant IDs -> we need at least 100, but wanted to have more data
# store each entry in dictionary
for rsid in variant_ids:
    rs_url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid}"
    rs_response = requests.get(rs_url)      # send GET request to NCBI

    rs_data = rs_response.json()    # store returned data

    # extract data
    try:
        placement = rs_data["primary_snapshot_data"]["placements_with_allele"][0]   # access first placement of variant (as they can have more than once placement on different assemblies)
        spdi = placement["alleles"][0]["allele"]["spdi"]    # spdi = Sequence Position Deletion Insertion (variant in NCBI)
        chrom = placement["seq_id"]     # choromosome ID
        pos = spdi["position"]      # genomic position
        ref = spdi["deleted_sequence"]      # reference allele
        alt = spdi["inserted_sequence"]     # alternate allele
    
        # determine variant type
        # comparing reference allele length to alternate allele length
        if len(ref) == 1 and len(alt) == 1:
            vtype = "SNP"
        elif len(ref) < len(alt):
            vtype = "Insertion"
        else:
            vtype = "Deletion"
        
        # store all variant info into a dictionary and append variant_records list from above
        variant_records.append({
            "rsid": rsid,
            "chromosome": chrom,
            "position": pos,
            "variant_type": vtype,
            "reference_allele": ref,
            "alternate_allele": alt
        })

    except (KeyError, IndexError, requests.RequestException):      # skip variants with missing fields or request issues
        continue

# write CSV
csv_file = "data.csv"      # new file (data)
with open(csv_file, "w", newline="") as f:
    # write CSV file with dictionaries as rows with our specified column headers
    writer = csv.DictWriter(f, fieldnames=[
         "rsid", "chromosome", "position", "variant_type", "reference_allele", "alternate_allele"                 
    ])
    writer.writeheader()    # write header row first
    for record in variant_records:      # write each variant dictionary as another row
        writer.writerow(record)

print(f"Variant data has been saved to {csv_file}") # print statement to show CSV has been made

# show preview of CSV as a dataframe
# load variant data from CSV into pandas dataframe
df = pd.read_csv("data.csv")

# print head of dataframe
print("\nPreview of variant data (pandas):")
print(df.head())
print(df.tail())



# STEP 2 - data processing & analysis
# Tasks:
    # Perform at least two different analyses on the dataset
        # Analyses:
        # 1) Variant Type Distribution Across APOE Gene (variant composition)
        # 2) Variant Densities Across APOE Gene (genomic location)
        # 3) Transition vs. Transversion Mutation Ratio In APOE Gene (mutation mechanism)
            # NOTE: see analysis_report.md for reasoning as to why this code was not included in this script

# ANALYSIS 1: Variant Type Distribution
    # What proportion of APOE variants are single nucleotide polymorphisms (SNPs) vs. insertions vs. deletions?

# ANALYSIS 1A: variant type counts with loops
# define function variant_type_distribution and pass parameters variant_records
    # count the number of SNPs, insertions, and deletions across the APOE gene
    # pass the list of dictionaries from step 1
    # return dictionary count of each variant type
def variant_type_distribution(variant_records):
    # set counts as 0 to hold counts for each variant type
    counts = {
        "SNP": 0,
        "Insertion": 0,
        "Deletion": 0
    }

    # run through variant_records and get values associated with variant_type key (SNP, insertion, deletion)
    for variant in variant_records:
        vtype = variant.get("variant_type")
        if vtype in counts:
            counts[vtype] += 1  # add one to count of whatever variant type the record is

    return counts

# run analysis 1
variant_type_counts = variant_type_distribution(variant_records)

# final print statements 
# count items in variant_type_counts and sort by vtype
print("\n\nStep 2: Analysis\nAnalysis 1A: Variant Type Distribution")
for vtype, count in variant_type_counts.items():
    print(f"{vtype}: {count}")

# ANALYSIS 1B: analysis with numpy
# pandas-based variant type distribution
variant_type_counts_pd = df["variant_type"].value_counts()

# variant type percentages with numpy
variant_type_percent = (
    variant_type_counts_pd / variant_type_counts_pd.sum()
) * 100

# print statements
print("\nAnalysis 1B: Variant Type Percentages")
print(variant_type_percent.round(2))


# ANALYSIS 2: Variant Densities 
    # Are variants evenly distributed across the APOE locus?

# ANALYSIS 2A: variant density with loops
# define function variant_density and pass it variants with bin_size=500
    # groups variants into genomic bins to analyze the positional density across the APOE gene
    # pass the list of dictionaries from step 1
    # return genomic bins dictionary 
        # NOTE: chose bin size of 500 after testing a few - bin size of 500 (in my case) provides four outputs, but bin size can be changed

def variant_density(variants, bin_size=500):
    # set bins dictionary to 0
    bins = {}

    # loop through each variant in the records list
    for v in variants:
        pos = v["position"]     # extract genomic position of the variant
        if pos == "":       # skip any variant with a missing position
            continue
        
        # find start position of the bin for said variant 
            # pos // bin_size tells how many full bins fit before this position
            # * bin_size gives start position of the bin for particular variant 
        bin_start = (pos // bin_size) * bin_size

        # create a readable label for the bin = start to end of 500bp bin 
        bin_label = f"{bin_start}-{bin_start + bin_size}"

        # count variants in the bin
        bins[bin_label] = bins.get(bin_label, 0) + 1

    return bins

# run analysis 2S
position_density = variant_density(variant_records)

# final print statements
print("\nAnalysis 2A: Variant Density by Genomic Region")
for region, count in position_density.items():
    print(f"{region}: {count}")

# ANALYSIS 2B: variant position stats with numpy
positions = df["position"].dropna()

# print stats
print("\nAnalysis 2B: Variant Position Statistics")
print("Mean position:", int(np.mean(positions)))
print("Minimum position:", int(np.min(positions)))
print("Maximum position:", int(np.max(positions)))

# ANALYSIS 2C: quantify clustering with pandas and numpy variance
# set bin size to match earlier
bin_size = 500

# make dataframe of bin counts
df["bin_start"] = (df["position"] // bin_size) * bin_size
bin_counts_pd = df["bin_start"].value_counts().sort_index()

# define function variant_clustering and pass it bin_counts (makes it reusable)
# calculate mean, variance, and standard deviation of variant density
def variant_clustering(bin_counts):
    mean = np.mean(bin_counts)
    variance = np.var(bin_counts)
    std_dev = np.std(bin_counts)

    return mean, variance, std_dev

# store values by calling variant_clustering
mean_density, density_variance, std_dev = variant_clustering(bin_counts_pd.values)

# print statements and values
print("\nAnalysis 2C: Quantify Variant Clustering")
print("Mean variants per bin:", mean_density)
print("Variance in variant density (variants per bin squared):", round(density_variance, 2))
print("Standard deviation of variant density (variants per bin):", round(std_dev, 2))



# STEP 3 - data visualization
# Tasks:
    # Use matplotlib or seaborn to generate at least two publication-quality visualizations representing data analysis results
        # Visualization 1: Pie Chart of Analysis 1
        # Visualization 2: Bar Chart of Analysis 2
    
# VISUALIZATION 1 - Variant Type Distribution Across the APOE Gene
# set labels and sizes for the pie chart (use info from variant_type_counts function)
labels = list(variant_type_counts.keys())
sizes = list(variant_type_counts.values())
color_map = {
    "SNP": "cornflowerblue",
    "Insertion": "steelblue",
    "Deletion": "indianred"
}
colors = [color_map[v] for v in labels]

plt.figure(figsize=(6,4))       # set figure size (longer for a bar chart)
bars = plt.bar(labels, sizes, color=colors)       # main bar chart
plt.xlabel("Variant Type")     # set axis labels
plt.ylabel("Number of Variants")      # set axis labels
plt.title("Variant Type Distribution Across the APOE Gene (500bp bins)")      # set title

# add labels above each bar
for bar in bars:
    height = bar.get_height()       # get bar height
    plt.text(
        bar.get_x() + bar.get_width() / 2,      # center the label horizontally above the bar
        height,
        str(height),
        ha="center",
        va="bottom"     # place text just above bar
    )

plt.tight_layout()
plt.savefig("figure1.png", dpi=300)       # save as a PNG
plt.close()

print("\n\nStep 3: Visualizations\nFigure has been saved as figure1.png")

# VISUALIZATION 2 - Variant Densities Across the APOE Gene
# sort bins by genomic position
sorted_bins = sorted(
    position_density.items(),       # turn dictionary into sequence of (key, value) pairs
    key=lambda x: int(x[0].split("-")[0])       # split bin label string into two parts, convert bin starting number (string) to a number, sort by these starting numbers
)

regions = [b[0] for b in sorted_bins]       # extract bin labels (bp-bp)
counts = [b[1] for b in sorted_bins]        # extract variant couns for each bin
colors = ["cornflowerblue", "darkorange", "seagreen", "indianred"]

plt.figure(figsize=(8,4))       # set figure size (longer for a bar chart)
bars = plt.bar(regions, counts, color=colors)       # main bar chart
plt.xlabel("Genomic Position (bp)")     # set axis labels
plt.ylabel("Number of Variants")      # set axis labels
plt.title("Variant Density Distribution Across the APOE Gene (500bp bins)")      # set title

# add labels above each bar
for bar in bars:
    height = bar.get_height()       # get bar height
    plt.text(
        bar.get_x() + bar.get_width() / 2,      # center the label horizontally above the bar
        height,
        str(height),
        ha="center",
        va="bottom"     # place text just above bar
    )

plt.tight_layout()
plt.savefig("figure2.png", dpi=300)       # save as a PNG
plt.close()

print("Figure has been saved as figure2.png")

# VISUALIZATION 3: Boxplot of Variant Density Per Bin
plt.figure(figsize=(4,4))
plt.boxplot(bin_counts_pd.values, vert=True)
plt.ylabel("Variants per 500 bp bin")
plt.title("Distribution of Variant Density\nAcross APOE Genomic Bins")

plt.tight_layout()
plt.savefig("figure3.png", dpi=400)     # save as a PNG
plt.close()

print("Figure has been saved as figure3.png\n")