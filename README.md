# CRISPR_targ
program to: 1. extract endogenous candidate Cas9/gRNA binding sites from a target sequence, and 2. screen them for orthogonality to context sequence (ex.: genome) 

version 1.0 (Aug. 21, 2015) 
-runs in MATLAB (built in R2014a)

-main script is CRISP_targ.m
-currently requires FASTA file(s) saved locally (input the full path to the file if it's not in same directory as the scripts). Future versions should be able to pull the sequences off genbank
-uses subroutines GG_library.m  (for 1.) and scan_SeqC.m (for 2.). 
