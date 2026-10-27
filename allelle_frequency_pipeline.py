#Allele Frequency Pipeline 
#This script was written to assist calculations with my dissertation at Florida State University: "Using Ecology and Evolutionary Biology to Explore the Diversity of Marine Nematodes"

#Author: H. Joseph Horacek, III
#Year: 2026
#This python script acts as a pipeline that takes diploid sequence data and calculates allele frequencies, the number of SNPs in each locus, and the allele counts per locus
#see the README file for more details

# Required package: pandas 
#
# This script reads an Excel file (.xlsx) containing four columns: Individual, Locus, Allele, Sequence 
# Individual --> The specimen that was sequenced
# Locus --> Genetic locus of the sequence
# Allele --> In this input file, this column identifies the two sequence copies observed for each diploid individual at a locus.  They are labeled 0 and 1. These are referred to as "alleles" in the sequencing/genotyping output.
#       #IMPORTANT: This is different from the "Allele" column in the allele-frequency output below. There, "Allele" refers to a unique sequence variant observed within a locus and is numbered 1, 2, 3...
# Sequence --> The DNA sequence for the Illumina allele
#
# An example input file is provided in the ExampleData folder.
#################

#Assumptions
#Sequences should already be aligned.
#Sequences should all be in the same orientation.
#Sequences should be the same length (per locus; different loci will have different sequence length)
#The script does not perform alignment.

import pandas as pd
from collections import Counter
import csv
#########################

# === CONFIG ===
input_excel = "loci-dummy.xlsx" ###this is the input file; see example data for layout of file

#Please note: for this to operate as the input file, the python file must be in the same folder as the input file
#   Also, the output files will be in the same folder 

#Input Excel file 
# Must contain at least four columns: Individual, Locus, Allele, Sequence 
#Individual	Locus	Allele	Sequence
#Specimen1	1	    0	    TGCAGATC....
#Specimen1	1	    1	    TGCAGATC....

#Please note that "allele" specifically refers to "Illumina alleles" which refer to the two sequence copies associated with an individual, representing the homologous chromosomes (diploid data).
#Alleles here are label either "0" or "1"; these are not the count of alleles but the label 
# Column names must be exactly:
#   Individual
#   Locus
#   Allele 
#   Sequence
#
# See the folder "ExampleData" in this repository for an example.

#Output1
output_csv = "allele_frequencies.csv" ###output for allele frequencies
## Allele frequencies for every locus
# Columns:
#   Locus 
#   Allele
#       [e.g. 1, 2, 3, etc.;these are unique sequence variants within the locus and are NOT the 0/1 Illumina allele labels]
#   Sequence
#   Count
#   Frequency

#Output2
allele_count_xlsx = "allele_counts_per_locus.xlsx" ###output for allele counts per locus
#Number of alleles observed per locus

#Output3
snp_count_xlsx = "snps_per_locus.xlsx" ###output for the number of SNPs per locus
#Number of SNP positions among unique alleles for each locus 


# === STEP 1: Load Excel ===
df = pd.read_excel(input_excel)

# === STEP 2: Calculate Allele Frequencies ===
output_rows = []
allele_counts_summary = []
snp_counts_summary = []

def count_snps(seqs):
    """Count SNP positions (differences across sequences)"""
    seqs = list(set(seqs))
    if len(seqs) < 2:
        return 0
    min_len = min(len(s) for s in seqs)
    snp_count = 0
    for i in range(min_len):
        bases = set(seq[i] for seq in seqs)
        if len(bases) > 1:
            snp_count += 1
    return snp_count

for locus, group in df.groupby("Locus"):
    allele_counts = Counter(group["Sequence"])
    total = sum(allele_counts.values())

    # Frequency table
    for i, (seq, count) in enumerate(
        sorted(allele_counts.items(), key=lambda x: x[1], reverse=True), 1
    ):
        freq = count / total
        output_rows.append({
            "Locus": locus,
            "Allele": i,
            "Sequence": seq,
            "Count": count,
            "Frequency": round(freq, 5)
        })

    # Allele count summary
    allele_counts_summary.append({
        "Locus": locus,
        "Num_Alleles": len(allele_counts)
    })

    # SNP count summary
        #Counts SNP positions among unique alleles.
        #Only unique sequences are compared.
        #If multiple bases occur at a nucleotide position, that position is counted as one SNP.
    snp_counts_summary.append({
        "Locus": locus,
        "Num_SNPs": count_snps(allele_counts.keys())
    })

# === STEP 3: Save Output Files ===

# Allele frequency table
with open(output_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Locus", "Allele", "Sequence", "Count", "Frequency"])
    writer.writeheader()
    writer.writerows(output_rows)

# Allele counts summary
pd.DataFrame(allele_counts_summary).to_excel(allele_count_xlsx, index=False)

# SNP counts summary
pd.DataFrame(snp_counts_summary).to_excel(snp_count_xlsx, index=False)

print(f"✅ Outputs written to:\n- {output_csv}\n- {allele_count_xlsx}\n- {snp_count_xlsx}")



##See README for more details regarding this script and methodology 

