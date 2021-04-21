# crc_subgroup
This repository contains scripts for querying the TCGA database and extracting data for COAD & READ projects. Files of patients from these projects having transcriptome data (HTseq.count) 
for "solid tissue normal" and "primary tumor" are downloaded and afterwards a differential geneexpression analysis will be performed in R with the DSEq2 pipline.

READ_final_V2_compact.py: Queries the READ project using python
COAD_final_V2_compact.py: Queries the COAD project using python
READ_htseqCounts_Analysis_V2.r: Analysis DGE of READ files between normal tissue and primary tumor using R
