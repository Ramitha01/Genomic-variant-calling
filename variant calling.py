#install GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip

#unzip the   GATK
unzip gatk-4.4.0.0.zip

#Navigate to the GATK directory:
cd gatk-4.4.0.0

(a)Variant Calling 
Add Read Groups (GATK):
gatk AddOrReplaceReadGroups \
  -I alignment/sample.sorted.bam \
  -O alignment/sample.RG.bam \
  -RGID sample \
  -RGLB LIB1 \
  -RGPL ILLUMINA \
  -RGPU UNIT1 \
  -RGSM sample

(b)Mark Duplicates:
gatk MarkDuplicates \
  -I alignment/sample.RG.bam \
  -O alignment/sample.dedup.bam \
  -M alignment/sample.metrics.txt

(c)Index BAM:
samtools index alignment/sample.dedup.bam

(d)Base Quality Score Recalibration (BQSR):
Create recalibration table:
gatk BaseRecalibrator \
  -I alignment/sample.dedup.bam \
  -R hg38.fa \
  --known-sites dbsnp_146.hg38.vcf.gz \
  --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O alignment/sample.recal.table

(e)Apply BQSR:
gatk ApplyBQSR \
  -R hg38.fa \
  -I alignment/sample.dedup.bam \
  --bqsr-recal-file alignment/sample.recal.table \
  -O alignment/sample.BQSR.bam

(f)Variant Calling (GATK HaplotypeCaller):
gatk HaplotypeCaller \
  -R hg38.fa \
  -I alignment/sample.BQSR.bam \
  -O variants/sample.raw.vcf.gz

#output:sample_raw.vcf.gz

Variant Filtering (Hard filters):
gatk VariantFiltration \
  -R hg38.fa \
  -V variants/sample.raw.vcf.gz \
  --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
  -O variants/sample.filtered.vcf.gz

#output:sample.filtered.vcf.gz



