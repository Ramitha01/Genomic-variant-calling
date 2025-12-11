Variant Calling and Functional Analysis Workflow

This repository documents a complete workflow for processing whole-genome or whole-exome sequencing data from raw FASTQ files to functionally interpreted variant sets. The workflow focuses on quality control, alignment, variant calling, annotation, and downstream biological interpretation.

The pipeline uses the following tools:
FastQC, Trim Galore, BWA-MEM, Samtools, GATK, ANNOVAR, and R.

1. Quality Assessment

The workflow begins by evaluating raw FASTQ files using standard QC metrics. This step helps identify issues such as low base quality, adapter contamination, uneven GC content, and overrepresented sequences. These metrics guide subsequent trimming and preprocessing decisions.

2. Adapter and Quality Trimming

Adapter sequences and low-quality regions are removed from the raw reads to improve alignment accuracy. Trimming ensures that only high-quality bases are used for downstream processes, reducing false alignments and improving variant calling reliability. Paired-end structure is maintained throughout.

3. Reference Preparation

The chosen reference genome (hg38) is prepared by generating all required index files. These include the BWA index for alignment, a FASTA index, and a sequence dictionary for GATK. Proper reference preparation ensures efficient and accurate read mapping and variant calling.

4. Read Alignment

Trimmed reads are aligned to the reference genome using the BWA-MEM algorithm, producing alignment files (BAM). These files contain the genomic positions of each read and serve as the foundation for all downstream GATK processing and variant discovery.

5. Post-Alignment Processing

Aligned reads are processed to meet best-practice standards. This includes sorting BAM files, ensuring correct ordering by genomic coordinate, and preparing the data for GATK-based variant calling. These steps improve data organization and enhance accuracy.

6. Variant Calling (GATK Workflow)

Variant calling in this pipeline follows GATK Best Practices and includes several preparatory steps required before generating reliable variant calls. All stages are performed using GATK.

a. Read Group Assignment

Read groups are added to the BAM files to embed metadata describing sample origin, library details, sequencing platform, and flowcell information. This metadata is essential for GATK tools to correctly interpret and process reads, especially when multiple samples or libraries are involved.

b. Duplicate Marking

PCR duplicates are identified and flagged to prevent them from influencing variant discovery. Duplicate reads do not represent independent biological observations and can artificially inflate read support for false variants. Marking duplicates reduces technical noise and improves the reliability of variant calls.

c. Base Quality Score Recalibration (BQSR)

BQSR corrects systematic biases in base quality scores using known variant databases. By recalibrating base qualities, sequencing machine–specific errors are reduced, resulting in more accurate modeling of true variants versus sequencing artifacts. This step significantly enhances variant-calling precision.

d. Variant Calling in GVCF Mode

Once all preprocessing steps are complete, variants are called using HaplotypeCaller in GVCF mode. This generates a per-sample genomic VCF (GVCF) that captures information across all genomic positions, including both variant and non-variant sites. The GVCF format supports accurate filtering and ensures a consistent representation of variants across samples.

7. Variant Annotation (ANNOVAR)

Annotated variant tables are generated using ANNOVAR. Multiple annotation sources are integrated, including gene models, clinical databases, functional prediction algorithms, population allele frequencies, cytogenetic mapping, and dbSNP identifiers. The resulting CSV annotation file provides detailed biological and functional context for each variant.

8. Tiered Variant Filtering in R

Annotated variants undergo a structured filtering approach to identify the most biologically relevant candidates.

a. Exonic Variant Extraction

Only variants located within protein-coding exons are retained for analysis, ensuring focus on regions most likely to affect protein function.

b. Functional Impact Filtering

Among exonic variants, only nonsynonymous SNVs are selected. These variants result in amino acid changes and are more likely to contribute to functional alterations and disease phenotypes.

c. Pathogenicity Filtering (SIFT)

Nonsynonymous SNVs are further filtered using SIFT scores. Only variants with SIFT ≤ 0.05, indicating predicted deleterious effect, are retained. This produces a refined set of high-confidence potentially damaging variants.

9. Gene List Preparation

Gene symbols corresponding to filtered variants are extracted and converted into Entrez Gene IDs. This ensures compatibility with enrichment analysis tools and standardized pathway databases.

10. Functional Enrichment Analysis

To interpret the biological significance of the filtered gene set, multiple enrichment analyses are performed.

GO Enrichment

Highlights enriched biological processes, cellular components, and molecular functions represented in the variant-affected genes.

KEGG Pathway Analysis

Identifies affected signaling pathways and metabolic routes with potential relevance to disease mechanisms.

Reactome Pathway Analysis

Provides mechanistic insight into molecular interactions and regulatory networks impacted by the variants.
All enrichment results are exported as spreadsheets and visualized using dot plots for ease of interpretation.

11. SNP–Gene Summary Table

A concise summary table is created linking each retained variant’s dbSNP identifier, gene symbol, and SIFT score. This serves as a practical summary for reporting, downstream investigation, or multi-omics integration.

12. Output Overview

The workflow generates:

1.	Quality-assessed and trimmed read files
2.	Aligned and fully processed BAM files
3.	Per-sample GVCF files generated by HaplotypeCaller
4.	Fully annotated variant tables
5.	Exonic, nonsynonymous, and SIFT-filtered variant lists
6.	GO/KEGG/Reactome enrichment results
7.	Pathway dot plots
