#Annotation:
#Convert VCF to ANNOVAR format:
perl annovar/convert2annovar.pl \
  -format vcf4 variants/sample.filtered.vcf.gz \
  -outfile annovar/sample.avinput \
  -includeinfo

#Download Required Databases
# Download RefSeq (Gene Names)
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

# Download Ensembl (Alternative Gene Names)
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/

# Download ClinVar (Disease Associations)
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20250721 humandb/

# Download dbSNP (Known SNPs)
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 humandb/

#Download cytoBand
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cytoBand humandb/

#Download 
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_exome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad211_genome humandb/

#Run ANNOVAR annotation:
perl table_annovar.pl VIZ0011.avinput humandb/ \
  -buildver hg38 \
  -out VIZ0011_annotated \
  -remove \
  -protocol refGene,ensGene,clinvar_20250721,dbnsfp47a,avsnp151,cytoband,gnomad211_exome \
  -operation g,g,f,f,f,f,f \
  -nastring . \
  -polish \
  -otherinfo \
  -csvout
