# IBDGem-LRs

This repository contains code used for the study *"On forensic likelihood ratios from low-coverage sequencing."* It is structured as follows:

## 1. Scripts

This directory contains all the scripts used to simulate the different scenarios mentioned in our paper.

   - **Diploid LD**: Contains scripts for simulating complete and perfect linkage disequilibrium scenarios in the diploid setting investigated in the paper.
   - **Diploid no LD**: Contains scripts for simulating diploid genotypes derived under neutral, CEU, or YRI site frequency spectra (SFS).
   - **Fraction of Unique Genotypes**: Contains scripts for working with the 1kGP VCF file to obtain results shown in Figure S5.
   - **Haploid**: Contains scripts for simulating haploid genotypes and running both complete and perfect LD scenarios.
   - **IBDGem Bash Scripts**: Contains the scripts used to run IBDGem on the simulated data.
   - **Recommendations**: Contains scripts for simulating scenarios in which IBDGem produces approximately correct outputs.
   - **Robustness to Phase**: Contains scripts to simulate two scenarios:
      1. Shuffled haplotypes in the reference panel.
      2. Inclusion of the target individual's haplotypes in the reference panel in either the same individual or across two different individuals.

## 2. Unique Genotype Stats

This directory contains the fractions of unique genotypes across the 100 different starting locations on chromosome 1, generated with different settings.

## 3. 1kGP

This directory contains output from the 1kGP analysis, including allele frequency distributions for each of the CEU and YRI subpopulations.

## Contact

For questions, issues, or contributions, please reach out through the GitHub repository's issue tracker or directly via email.
