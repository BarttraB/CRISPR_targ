# CRISPR_targ

Program to: 1. extract endogenous candidate Cas9/gRNA binding sites from a target sequence, and 2. screen them for orthogonality to context sequence (ex.: genome) 

*Matlab branch:*

Made in R2014a; Version 1.3 (Aug. 24, 2015) 

-main script to run is CRISP_targ

-uses subroutines GG_library (for 1.) and scan_SeqC (for 2.). 

-currently requires FASTA file(s) saved locally (input the full path to the file if it's not in same directory as the scripts). Future versions should be able to pull the sequences off genbank

