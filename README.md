This repository contains step-by-step instructions, code and data for:  
Maddie E. James, Travis G. Britton, Jonathan D. Mitchell, Ben Halliwell, Barbara Holland, Ian J. Wright, Daniel Ortiz-Barrientos. Pervasive linked selection and non-convergent aridity adaptation in *Eucalyptus*.

# Genetics analysis

## Bioinformatics

Leaf samples were sent to the Australian Genome Research Facility for DNA extractions and library preparation with the Illumina DNA Prep M Library Prep kit. Samples were pooled in equimolar concentrations and sequenced on two lanes of the NovaSeq 600 S4, 300 cycle with 150bp paired-end sequencing.

### Quality filtering and trimming
For each sample, raw reads were quality-trimmed with `fastp` v0.23.4, using: minimum base quality 10, maximum uncalled bases 50, and minimum read length 50 bp. An HTML report was generated for each sample.

```
fastp --in1 forward.fastq.gz --in2 reverse.fastq.gz \
--out1 trimmed_R1.fastq.gz --out2 trimmed_R2.fastq.gz \
-q 10 -u 50 -l 50 -w 1 -h sample_fastp.html
```

### Read alignment
Reads were mapped to the *Eucalyptus grandis* v2.0 reference genome (https://phytozome-next.jgi.doe.gov/info/Egrandis_v2_0) using `BWA-MEM` v0.7.19. Read groups were added during mapping. The SAM output was piped into `SAMtools` v1.19 for coordinate sorting and BAM creation.

```
bwa mem -t 12 -M -R "@RG\tSM:sample\tID:sampleID\tLB:lib1\tPL:ILLUMINA\tPU:barcode" reference.fasta \
trimmed_R1.fastq.gz trimmed_R2.fastq.gz | \
samtools sort -@ 12 -o sample_srt.bam
```

### Alignment statistics
Alignment statistics were generated with `SAMtools` v1.19.

```
samtools flagstat sample_srt.bam > sample_srt_stats.txt
```

### Cleaning BAM files
`Picard` v3.1.1 CleanSam was used to clean each BAM (to set soft-clipping for reads beyond end of the reference, and MAPQ to 0 for unmapped reads).

```
java -Xmx8g -jar picard.jar CleanSam \
INPUT=sample_srt.bam \
OUTPUT=sample_cln.bam
```

### PCR duplicate marking and indexing
PCR duplicates were marked with `Picard` MarkDuplicates. The resulting BAM files were indexed with `SAMtools`, and per-sample coverage was calculated with SAMtools coverage.

```
java -Xmx4g -jar picard.jar MarkDuplicates \
INPUT=sample_cln.bam \
OUTPUT=sample_PCRm_cln.bam \
METRICS_FILE=sample_PCRm_cln_metrics

samtools index sample_PCRm_cln.bam

samtools coverage sample_PCRm_cln.bam --output sample_coverage.txt
```

### Variant calling
Per-sample variant discovery was performed with `GATK` v4.4.0.0 HaplotypeCaller in GVCF mode, retaining both variant and invariant sites. The example below shows one chromosome. Scripts were run in parallel across chromosomes and samples.

```
gatk --java-options "-Xmx8g" HaplotypeCaller \
--input sample_PCRm_cln.bam \
--output sample_Chr01_vrnt.g.vcf.gz \
--reference reference.fasta \
--emit-ref-confidence GVCF \
--intervals Chr01
```

Chromosome-level GVCFs for each sample were concatenated and indexed with `bcftools`.

```
bcftools concat --threads 6 \
sample_Chr01_vrnt.g.vcf.gz \
sample_Chrn_vrnt.g.vcf.gz \
-o sample_vrnt.g.vcf.gz
bcftools index --threads 6 -t sample_vrnt.g.vcf.gz
```

Individuals were combined into cohorts (i.e. one file per species) with `GATK` CombineGVCFs.

```
gatk CombineGVCFs \
-R reference.fasta \
--variant sample01_vrnt.g.vcf.gz \
--variant samplen_vrnt.g.vcf.gz \
-O species1.g.vcf.gz
```

### Joint genotyping 
For each species, joint genotyping was performed per chromosome with `GATK` GenotypeGVCFs, retaining variant and invariant sites. The example below shows one chromosome. Scripts were run in parallel across chromosomes and samples.

```
gatk --java-options "-Xmx8g" GenotypeGVCFs \
-R reference.fasta \
-V species1.g.vcf.gz \
-O species1_Chr01_ species1.g.vcf.gz \
--include-non-variant-sites \
--intervals Chr01
```

Chromosome-level cohort VCFs were concatenated and indexed with `bcftools`.

```
bcftools concat --threads 6 \
species1_Chr01_ species1.g.vcf.gz \
species1_Chrn_ species1.g.vcf.gz \
-o species1_allchrom.g.vcf.gz
bcftools index --threads 6 -t species1_allchrom.g.vcf.gz
```

All populations were merged with `bcftools` into a single file, then bgzip-indexed.

```
bcftools merge --threads 96 \
--file-list bcftoolsmerge.txt \
--output allPops.g.vcf.gz
bcftools index --threads 96 -t allPops.g.vcf.gz
```

Missing genotypes encoded as reference with DP=0 were corrected using `bcftools` +setGT. This was because GATK v4.4.0.0 encodes missing genotypes as 0/0 with a DP=0 in the format field for single-sample SNP calling with HaplotypeCaller, which can cause downstream issues. The corrected VCF was bgzip-compressed and indexed with `bcftools`.

```
bcftools +setGT allPops.g.vcf.gz -- -t q -n . -e 'FMT/DP=0' > allPops-Cor.g.vcf
bcftools view allPops-Cor.g.vcf -Oz -o allPops-Cor.g.vcf.gz
bcftools index -t allPops-Cor.g.vcf.gz
```

### SNP filtering
Site-level hard filters were applied with `GATK` VariantFiltration.

```
gatk VariantFiltration \
-V allPops-Cor.g.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 30.0" --filter-name "MQ30" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O allPops-CorVFnoQual.g.vcf.gz
```

Using `VCFtools` v0.1.17, sites failing GATK filters were removed, and variants were retained if genotyped in at least 50% of individuals.

```
vcftools --gzvcf allPops-CorVFnoQual.g.vcf.gz \
--max-missing 0.5 \
--remove-filtered-all \
--recode --recode-INFO-all --stdout | gzip -c > allPops-CorVFnoQualMM.g.vcf.gz
```

Minimum depth per sample of 3 reads was applied; genotypes with fewer than 3 reads were recoded as missing.

```
vcftools --gzvcf allPops-CorVFnoQualMM.g.vcf.gz \
--minDP 3 \
--recode --recode-INFO-all --stdout | gzip -c > allPops-CorVFnoQualMMdp3.g.vcf.gz
```

To identify low-quality individuals (high missingness), per-individual missing data proportions were computed. The output file (out.imiss) has the fifth column as the proportion missing.

```
vcftools --gzvcf allPops-CorVFnoQualMMdp3.g.vcf.gz \
--missing-indv \
--out allPops-CorVFnoQualMMdp3MissingInds
```

<img src="images/Missing_data_per_individual.png" width="700" />

As all individuals contained less than 40% missing data, they were all retained.

Mean depth per site was calculated.

```
vcftools --gzvcf allPops-CorVFnoQualMMdp3.g.vcf.gz \
--site-mean-depth \
--out allPops-CorVFnoQualMMdp3MeanDepth.txt
```

Sites were retained if mean depth was between 10 and 100, and only biallelic sites were kept.

```
vcftools --gzvcf allPops-CorVFnoQualMMdp3.g.vcf.gz \
--max-meanDP 100 \
--min-meanDP 10 \
--max-alleles 2 \
--recode --recode-INFO-all \
--out allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz
```

To ensure representation across populations, we created a 'keep' list, where individuals were grouped by species using a two-column file (individual_name, species_id). The proportion of missing data per site was computed for each species (produces *.lmiss files; column 6 is proportion missing). 

```
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep dealba_keep   --missing-site --out dealba
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep largiflo_keep --missing-site --out largiflo
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep moluc_keep    --missing-site --out moluc
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep nitida_keep   --missing-site --out nitida
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep regnans_keep  --missing-site --out regnans
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep siderox_keep  --missing-site --out siderox
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep stricta_keep  --missing-site --out stricta
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep vicina_keep   --missing-site --out vicina
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf --keep viridis_keep  --missing-site --out viridis
```

A bad loci list was built per species (>70% missing in any species), keeping CHR and POS columns. 

```
cat dealba.lmiss largiflo.lmiss moluc.lmiss nitida.lmiss regnans.lmiss siderox.lmiss stricta.lmiss vicina.lmiss viridis.lmiss \
| awk 'NR>1 && $6 > 0.7 {print $1"\t"$2}' > badloci0.7
```

These positions were excluded from the VCFfile.

```
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10biallele.g.vcf.gz.recode.vcf \
--exclude-positions badloci0.7 \
--recode --recode-INFO-all \
--out allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70
```

Sites were removed with >20% overall missing data.

```
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70.recode.vcf \
--max-missing 0.8 \
--recode --recode-INFO-all \
--out allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70MM80
```

Indels were removed.

```
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70MM80.recode.vcf \
--remove-indels \
--recode --recode-INFO-all \
--out allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70MM80noIndelsVCFtools
```

Final dataset (allPops.vcf): 160,686,818 variant and invariant sites; 12,608,766 variant sites with MAF ≥ 0.05. 

## Phylogeny

A coalescent-based species tree was inferred across all species. 

### Gene alignments
Per-gene VCFs were converted to PHYLIP alignments with `vcf2phylip` v2.0 (IUPAC codes for heterozygous sites).

```
python3 vcf2phylip.py -i allPopsGene${i}.vcf -p
```

### Gene trees
One tree per gene (all individuals) was inferred with `IQ-TREE` v2.3.6, run once per gene (e.g. with GNU `parallel`).

```
iqtree2 -s allPopsGene${i}.phy -m GTR+I+G -nt 1 -czb --redo
```

### Collapsing near-zero branches
Branches shorter than 0.000001 were collapsed to polytomies using `collapseEdges` from the R package `MSCquartets`.

```r
library(ape)
library(MSCquartets)

trees <- read.tree("GeneTrees.txt")
trees_collapsed <- lapply(trees, function(tr) collapseEdges(tr, 1e-6))
class(trees_collapsed) <- "multiPhylo"
write.tree(trees_collapsed, file = "GeneTrees_collapsed.txt")
```

### Species tree
The species tree was inferred from the collapsed gene trees with `ASTRAL-IV`.

```
astral4 -t 25 -i GeneTrees_collapsed.txt -o ASTRALSpeciesTree.txt
```

## Population structure

Principal components analysis (PCA) on SNPs with MAF ≥ 0.05, LD-pruned in `PLINK` v1.9 (window of 50 SNPs, step of 10 SNPs, r² < 0.1), giving 428,330 SNPs.

### LD pruning and PCA
```
plink --vcf allPopsMAF0.05NoScaffolds.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out allPopsUnlinked

plink --vcf allPopsMAF0.05NoScaffolds.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --extract allPopsUnlinked.prune.in --make-bed --pca --out allPopsUnlinked.PCA
```

## Genome-wide landscapes of diversity and divergence (20kb windows)

Genome-wide FST and *d*XY (between species) and π and Tajima's D (within species) were calculated for the three species pairs with 10 individuals per species: the two wet-dry contrasts of *E. moluccana*–*E. largiflorens* and *E. dealbata*–*E. vicina*, and the wetter-wet contrast of *E. nitida*–*E. regnans*. FST, *d*XY and π were calculated with 20kb non-overlapping sliding windows in `pixy` v1.2.10, which accounts for monomorphic sites when calculating diversity and divergence statistics. 

```
pixy --stats pi fst dxy \
--vcf allPops.vcf.gz \
--populations allPops_popfile.txt \
--output_prefix 20000 \
--window_size 20000 \
--n_cores 6 \
```

Tajima’s D was calculated for each population in VCFtools using the same 20kb non-overlapping sliding windows. 

```
vcftools --vcf allPops.vcf --keep dealba.txt --TajimaD 20000 --out dealba_tajimasD_20000
vcftools --vcf allPops.vcf --keep largiflo.txt --TajimaD 20000 --out largiflo_tajimasD_20000
vcftools --vcf allPops.vcf --keep moluc.txt --TajimaD 20000 --out moluc_tajimasD_20000
vcftools --vcf allPops.vcf --keep nitida.txt --TajimaD 20000 --out nitida_tajimasD_20000
vcftools --vcf allPops.vcf --keep regnans.txt --TajimaD 20000 --out regnans_tajimasD_20000
vcftools --vcf allPops.vcf --keep vicina.txt --TajimaD 20000 --out vicina_tajimasD_20000
```

We used custom R code [CombiningFiles-allThreePairs.R](Data%20files%20and%20code/CombiningFiles-allThreePairs.R) to combine the above population genetics statistics into a single file per population pair. 

Files are as follows:  
*E. moluccana*–*E. largiflorens*: [moluc_largiflo_merged.csv](Data%20files%20and%20code/moluc_largiflo_merged.csv)  
*E. dealbata*–*E. vicina*: [dealba_vicina_merged.csv](Data%20files%20and%20code/dealba_vicina_merged.csv)  
*E. nitida*–*E. regnans*: [nitida_regnans_merged.csv](Data%20files%20and%20code/nitida_regnans_merged.csv)

### Plotting landscapes
Genome-wide landscapes of FST, *d*XY, π and Tajima's D for each species pair. For visualisation only, *d*XY, π and Tajima's D were smoothed by taking the mean across 1.5Mb bins.

R code for plotting genome-wide landscapes of diversity and divergence for each species pair: [Plotting-Lanscapes.R](Data%20files%20and%20code/Plotting-Lanscapes.R).  

### Correlations of landscapes among species (20kb windows)
Pearson's *r* and Spearman's ρ correlations of π, Tajima's D, *d*XY and FST landscapes across non-overlapping 20kb windows. Code: [Correlation-20kb-Landscapes.R](Data%20files%20and%20code/Correlation-20kb-Landscapes.R).  

## Linked selection throughout the genome (1Mb windows)

To test for linked selection, π, FST and *d*XY were calculated in non-overlapping 1Mb windows for each species pair, alongside genome features in the same 1Mb windows.

### Population genetics statistics and genome features (1Mb windows)
π, FST and *d*XY in 1Mb windows with `pixy`(as above, with 1Mb windows). GC content from the *E. grandis* v2.0 reference was calculated with `bedtools` (as outlined below, with 1Mb windows). Gene density from the reference annotation using a custom R script (as outlined below); a gene was assigned to a window if more than 50% of its length overlapped the window. Recombination rate (cM/Mb) from the linkage map used to construct the *E. grandis* reference genome (Version 2 markers, GRA parent; Bartholomé et al. 2015).

### Residualising recombination rate
Recombination rate (natural log) was residualised against gene density and GC content using a generalised additive model (thin-plate smooth terms; `mgcv` in R). The residuals were used to assess relationships with mean π, FST and *d*XY.

```r library(mgcv)
data$log_rec <- log(data$rec_cM_Mb)
gam_model <- gam(log_rec ~ s(gene_count) + s(pct_gc), data = data, method = "REML") resid_rec <- residuals(gam_model) ```

### Correlations and plots (1Mb windows)
Pearson's *r* and Spearman's ρ across 1Mb windows for: FST vs mean π, *d*XY vs mean π, mean π vs gene density, mean π vs recombination rate, and relationships among genome features. This used the same code as for the 20kb windows (above). 

### Robustness check: outgroup normalisation
To control for window-to-window mutation rate variation, π was divided by each species' *d*XY to the outgroup *E. stricta*, and *d*XY for each pair was divided by the mean *d*XY of its two species to *E. stricta*. The analyses above were repeated with the normalised statistics. 

Data file: [allSpecies-1MbWindows.csv](Data%20files%20and%20code/allSpecies-1MbWindows.csv)

## Detection of genomic outliers (20kb windows)

Outliers were detected for each species pair using 20kb windows.

### Residual-based approach (top 1%)
*d*XY was regressed on mean π (linear model) and FST on mean π (LOESS, span = 0.6), and residuals extracted per window. Outlier windows were those in the top 1% for both residual FST and residual *d*XY. 

```r
dxy_residual <- residuals(lm(dxy ~ mean_pi))
fst_residual <- residuals(loess(fst ~ mean_pi, span = 0.6))
```

Files are as follows:  
*E. moluccana*–*E. largiflorens*: [moluc_largiflo_merged_outliers.csv](Data%20files%20and%20code/moluc_largiflo_merged.csv)  
*E. dealbata*–*E. vicina*: [dealba_vicina_merged_outliers.csv](Data%20files%20and%20code/dealba_vicina_merged.csv)  
*E. nitida*–*E. regnans*: [nitida_regnans_merged_outliers.csv](Data%20files%20and%20code/nitida_regnans_merged.csv)

### Sharing of outlier windows among pairs
Shared outlier windows between each pair of species pairs were tested with a circular-shift permutation (positions rotated by a random offset within each chromosome; 10,000 replicates). Aridity-specific sharing was summarised as Δ = E₁₂ − (E₁₃ + E₂₃)/2, with the null SD of Δ from circular-shift replicates. Code: [CircularShiftNull-OutlierSharing.R](Data%20files%20and%20code/CircularShiftNull-OutlierSharing.R).

### Assigning genes to outlier windows
A gene was considered to be within an outlier window if at least 25% of its length resided in the window. *Arabidopsis* orthologues were obtained from the annotated *E. grandis* reference genome. Code: [ExtractingGenes-OutlierWindows.R](Data%20files%20and%20code/ExtractingGenes-OutlierWindows.R).

## Association of outliers with genome features (20kb windows)

GC content for 20kb windows was obtained along the *E. grandis* reference genome with `bedtools`.

```
bedtools makewindows -g Eucalyptus_grandis.main_genome.scaffolds.fasta.fai -w 20000 > Eucalyptus_grandis.main_genome.scaffolds_windows_20kb.bed
bedtools nuc -fi Eucalyptus_grandis.main_genome.scaffolds.fasta -bed Eucalyptus_grandis.main_genome.scaffolds_windows_20kb.bed > Eucalyptus_grandis.main_genome.scaffolds_windows_20kb_GCcontent.txt
```

For each species pair, outlier windows were compared with background windows for π, recombination rate, gene density (genes per 20kb), GC content, FST and *d*XY. The test statistic was the median of the outlier windows, compared against a null of 10,000 random samples of background windows matched in number (two-sided empirical P-value). Code: [Outlier-Features.R](Data%20files%20and%20code/Outlier-Features.R)

## Genomic islands of differentiation (20kb windows)

Broad islands of differentiation (extended contiguous regions of elevated FST rather than individual outlier windows) were identified using the 20kb window FST values from `pixy`.

### Identifying islands
Islands were regions in the top 5% of the genome-wide FST distribution spanning at least ten consecutive 20kb windows, allowing up to two consecutive non-outlier windows within an island. Windows lacking data were excluded. Islands separated by only a few windows were merged (23 islands in total).

### Island features compared with background
Each island was compared with background windows (all windows not in any island) for π in each population, *d*XY, recombination rate, gene density and GC content using two-sided Wilcoxon rank-sum tests, with Benjamini-Hochberg correction within each feature and species pair. Code: [Islands-vs-Background.R](Data%20files%20and%20code/Islands-vs-Background.R)
