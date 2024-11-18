import pandas as pd
import numpy as np 

ceu_allele_frequencies = pd.read_csv('ceu_allele_frequencies.txt', sep='\s+', usecols=[4], header = None)
ceu_allele_freqs = ceu_allele_frequencies.iloc[:, 0]
frequency_counts = ceu_allele_freqs.value_counts().sort_index()
total_snps = frequency_counts.sum()
ceu_sfs = (frequency_counts / total_snps).values
np.savetxt('ceu_sfs.txt', ceu_sfs)

yri_allele_frequencies = pd.read_csv('yri_allele_frequencies.txt', sep='\s+', usecols=[4], header = None)
yri_allele_freqs = yri_allele_frequencies.iloc[:, 0]
frequency_counts = yri_allele_freqs.value_counts().sort_index()
total_snps = frequency_counts.sum()
yri_sfs = (frequency_counts / total_snps).values
np.savetxt('yri_sfs.txt', yri_sfs)
