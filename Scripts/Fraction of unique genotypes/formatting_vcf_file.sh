#!/bin/bash

# Filter to keep only biallelic sites
bcftools view -m2 -M2 -v snps output.vcf.gz -o biallelic_snps.vcf.gz

# Apply minor allele frequency filter on the biallelic sites
bcftools view -q 0.01:minor biallelic_snps.vcf.gz -o filtered_biallelic_snps.vcf.gz

# Query both filtered and unfiltered vcfs to create genotype matrices
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' biallelic_snps.vcf.gz > biallelic_genotype_matrix.txt &
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' filtered_biallelic_snps.vcf.gz > filtered_genotype_matrix.txt &
