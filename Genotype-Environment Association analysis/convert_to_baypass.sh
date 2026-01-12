
### This script converts allele frequencies from Plink format (per population) into BayPass format. Note that the ---keep-allele-order is required in Plink to keep ref/alt alleles consistent across the files (column A1 is ref, column A2 is alt).

# Define populations in order
populations=(dealba largiflo moluc nitida regnans siderox stricta vicina viridis)

echo "Converting PLINK frequency files to BayPass format (corrected REF/ALT)..."

# Get total SNPs from first file
total_snps=$(awk 'END {print NR-1}' dealba_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq)
echo "Processing $total_snps SNPs across ${#populations[@]} populations"

# Use awk to process all files simultaneously
awk '
BEGIN {
    # Population file names
    pop[1] = "dealba_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[2] = "largiflo_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq" 
    pop[3] = "moluc_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[4] = "nitida_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[5] = "regnans_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[6] = "siderox_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[7] = "stricta_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[8] = "vicina_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    pop[9] = "viridis_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq"
    npops = 9
}

# Skip header lines (FNR==1 for each file)
FNR==1 {next}

# Process data lines
{
    # Store data for each population
    snp_id = $2
    freq = $5      # This is MAF = frequency of A2 = frequency of ALT allele
    nchrobs = $6
    
    # Calculate allele counts
    # PLINK: A1=REF, A2=ALT, MAF=freq(A2)=freq(ALT)
    if (nchrobs == 0) {
        ref_count = 0
        alt_count = 0
    } else {
        alt_count = int(freq * nchrobs + 0.5)      # MAF = A2 frequency = ALT frequency
        ref_count = nchrobs - alt_count            # REF count is remainder
    }
    
    # Store counts for this population
    for (i = 1; i <= npops; i++) {
        if (FILENAME == pop[i]) {
            ref_counts[FNR, i] = ref_count  # REF count (A1)
            alt_counts[FNR, i] = alt_count  # ALT count (A2)
            break
        }
    }
}

END {
    # Output the matrix
    for (snp_line = 2; snp_line <= FNR; snp_line++) {
        row = ""
        for (pop_i = 1; pop_i <= npops; pop_i++) {
            row = row " " ref_counts[snp_line, pop_i] " " alt_counts[snp_line, pop_i]
        }
        print substr(row, 2)  # Remove leading space
    }
}
' "${populations[@]/%/_allelefreq_LD-50-5-0.2_chr1-11_top10FstRemoved.frq}" > baypass_allele_counts_MAF0.01_LD-50-5-0.2_chr1-11_top10FstRemoved.txt

echo "BayPass allele count matrix created: baypass_allele_counts_MAF0.01_LD-50-5-0.2_chr1-11_top10FstRemoved.txt"
echo "Matrix dimensions: $total_snps SNPs × $((${#populations[@]} * 2)) columns"

# Show first few lines for verification
echo -e "\nFirst 3 lines of output:"
head -3 baypass_allele_counts_MAF0.01_LD-50-5-0.2_chr1-11_top10FstRemoved.txt

# Verify column count
echo -e "\nColumn count check:"
head -1 baypass_allele_counts_MAF0.01_LD-50-5-0.2_chr1-11_top10FstRemoved.txt | wc -w

echo -e "\nConversion complete!"
echo "Population order: ${populations[*]}"
echo "Column pairs: REF ALT REF ALT ... (for each population in order)"
