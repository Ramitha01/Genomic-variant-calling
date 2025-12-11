#Set Up Environment
# Source the .bashrc file
source ~/.bashrc

# Create the working environment for genomics obesity analysis
conda create -n genomics_pipeline
conda activate genomics_pipeline

# Create project directory
mkdir -p genomics_pipeline/{raw,trimmed,fastqc,alignment,variants,annovar}

##FASTQC
# Run FastQC to check the quality of the raw data
fastqc -o /path/to/output /path/to/input.fastq

##TRIMMING
trim_galore --paired raw/input_1.fastq.gz raw/input_2.fastq.gz \
  --output_dir trimmed/

#Download the reference genome
wget ftp://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#unzip the reference genome 
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

##Reference Genome Indexing with BWA
bwa index reference_genome.fa

##Alignment with BWA MEM
bwa mem reference_genome.fa input_trimmed.fastq > aligned.sam

##Convert SAM to BAM with Samtools
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam -o sorted_aligned.bam
samtools index sorted_aligned.bam
