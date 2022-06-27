Oxford-CTS cohort analysis

Genomic data quality control described in 'Oxford-CTS Genomic QC.ipynb'
This pipeline was used to generate the genotype calls for rs62175241
Using available paired RNA-seq data (raw counts in file 'Oxford_CTS_raw_counts.csv'),  we performed eQTL analysis using the R script ('Oxford-CTS eQTL.R')

This script requires:

- Matrix mapping RNA-seq to genotyping unique IDs ('Translating sample to CTS_ID.txt'), also comprising age and gender for participants
- Genotype calls at rs62175241 SNP ('cts_dirc3_snp.csv')
- Principal components computed in this cohort using pruned SNPs ('cts_pcs.tsv'), derivation described in 'Oxford-CTS Genomic QC.ipynb'