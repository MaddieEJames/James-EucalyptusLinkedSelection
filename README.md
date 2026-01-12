# James-et-al.2026-Eucalyptus

This repository contains step-by-step instructions, code and data for:
Maddie E. James, Jonathan D. Mitchell, Travis Britton, Ben Halliwell, Barbara Holland, Ian Wright, Daniel Ortiz-Barrientos.
Genotypic and phenotypic convergence due to climatic adaptation in *Eucalyptus*.

# Trait analysis
[to complete]

# Genetics analysis

## Bioinformatics

Leaf samples were sent to the Australian Genome Research Facility for DNA extractions and library preparation with the Illumina DNA Prep M Library Prep kit. Samples were pooled in equimolar concentrations and sequenced on two lanes of the NovaSeq 600 S4, 300 cycle with 150bp paired-end sequencing.

### Quality filtering and trimming
For each sample, raw reads were quality-trimmed with ```fastp``` v0.23.4, using: minimum base quality 10, maximum uncalled bases 50, and minimum read length 50 bp. An HTML report was generated for each sample.

```
fastp --in1 forward.fastq.gz --in2 reverse.fastq.gz \
--out1 trimmed_R1.fastq.gz --out2 trimmed_R2.fastq.gz \
-q 10 -u 50 -l 50 -w 1 -h sample_fastp.html
```

### Read alignment
Reads were mapped to the *Eucalyptus grandis* v2.0 reference genome (https://phytozome-next.jgi.doe.gov/info/Egrandis_v2_0) using ```BWA-MEM``` v0.7.19. Read groups were added during mapping. The SAM output was piped into ```SAMtools``` v1.19 for coordinate sorting and BAM creation.

```
bwa mem -t 12 -M -R "@RG\tSM:sample\tID:sampleID\tLB:lib1\tPL:ILLUMINA\tPU:barcode" reference.fasta \
trimmed_R1.fastq.gz trimmed_R2.fastq.gz | \
samtools sort -@ 12 -o sample_srt.bam
```

### Alignment statistics
Alignment statistics were generated with ```SAMtools``` v1.19.

```
samtools flagstat sample_srt.bam > sample_srt_stats.txt
```

### Cleaning BAM files
```Picard``` v3.1.1 CleanSam was used to clean each BAM (to set soft-clipping for reads beyond end of the reference, and MAPQ to 0 for unmapped reads).

```
java -Xmx8g -jar picard.jar CleanSam \
INPUT=sample_srt.bam \
OUTPUT=sample_cln.bam
```

### PCR duplicate marking and indexing
PCR duplicates were marked with ```Picard``` MarkDuplicates. The resulting BAM files were indexed with ```SAMtools```, and per-sample coverage was calculated with SAMtools coverage.

```
java -Xmx4g -jar picard.jar MarkDuplicates \
INPUT=sample_cln.bam \
OUTPUT=sample_PCRm_cln.bam \
METRICS_FILE=sample_PCRm_cln_metrics

samtools index sample_PCRm_cln.bam

samtools coverage sample_PCRm_cln.bam --output sample_coverage.txt
```

### Variant calling
Per-sample variant discovery was performed with ```GATK``` v4.4.0.0 HaplotypeCaller in GVCF mode, retaining both variant and invariant sites. The example below shows one chromosome. Scripts were run in parallel across chromosomes and samples.

```
gatk --java-options "-Xmx8g" HaplotypeCaller \
--input sample_PCRm_cln.bam \
--output sample_Chr01_vrnt.g.vcf.gz \
--reference reference.fasta \
--emit-ref-confidence GVCF \
--intervals Chr01
```

Chromosome-level GVCFs for each sample were concatenated and indexed with ```bcftools```.

```
bcftools concat --threads 6 \
sample_Chr01_vrnt.g.vcf.gz \
sample_Chrn_vrnt.g.vcf.gz \
-o sample_vrnt.g.vcf.gz
bcftools index --threads 6 -t sample_vrnt.g.vcf.gz
```

Individuals were combined into cohorts (i.e. one file per species) with ```GATK``` CombineGVCFs.

```
gatk CombineGVCFs \
-R reference.fasta \
--variant sample01_vrnt.g.vcf.gz \
--variant samplen_vrnt.g.vcf.gz \
-O species1.g.vcf.gz
```

### Joint genotyping 
For each species, joint genotyping was performed per chromosome with ```GATK``` GenotypeGVCFs, retaining variant and invariant sites. The example below shows one chromosome. Scripts were run in parallel across chromosomes and samples.

```
gatk --java-options "-Xmx8g" GenotypeGVCFs \
-R reference.fasta \
-V species1.g.vcf.gz \
-O species1_Chr01_ species1.g.vcf.gz \
--include-non-variant-sites \
--intervals Chr01
```

Chromosome-level cohort VCFs were concatenated and indexed with ```bcftools```.

```
bcftools concat --threads 6 \
species1_Chr01_ species1.g.vcf.gz \
species1_Chrn_ species1.g.vcf.gz \
-o species1_allchrom.g.vcf.gz
bcftools index --threads 6 -t species1_allchrom.g.vcf.gz
```

All populations were merged with ```bcftools``` into a single file, then bgzip-indexed.

```
bcftools merge --threads 96 \
--file-list bcftoolsmerge.txt \
--output allPops.g.vcf.gz
bcftools index --threads 96 -t allPops.g.vcf.gz
```

Missing genotypes encoded as reference with DP=0 were corrected using ```bcftools``` +setGT. This was because GATK v4.4.0.0 encodes missing genotypes as 0/0 with a DP=0 in the format field for single-sample SNP calling with HaplotypeCaller, which can cause downstream issues. The corrected VCF was bgzip-compressed and indexed with ```bcftools```.

```
bcftools +setGT allPops.g.vcf.gz -- -t q -n . -e 'FMT/DP=0' > allPops-Cor.g.vcf
bcftools view allPops-Cor.g.vcf -Oz -o allPops-Cor.g.vcf.gz
bcftools index -t allPops-Cor.g.vcf.gz
```

### SNP filtering
Site-level hard filters were applied with ```GATK``` VariantFiltration.

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

Using ```VCFtools``` v0.1.17, sites failing GATK filters were removed, and variants were retained if genotyped in at least 50% of individuals.

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

[put missing distribution graph here!!! ]

As all individuals contained > XXX, then they were all retained. 





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

Per-species missingness was computed, and sites with >70% missing in any species were removed. Overall missingness >20% was also filtered out. Indels were removed.

```
vcftools --vcf allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70MM80.recode.vcf \
--remove-indels \
--recode --recode-INFO-all \
--out allPops-CorVFnoQualMMdp3MaxMD100MinMD10bialleleMpp70MM80noIndelsVCFtools
```

Final dataset: 160,686,818 variant and invariant sites; 12,608,766 variant sites with MAF ≥ 0.05.




## Convergence divergence models
[to complete]

## Population genetics analyses
[to complete]

## Genotype-Environment Association (GEA) analysis
[to complete]
