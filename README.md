# TRaCE
Transcript Ranking and Canonical Election
## running TRaCE
```./TRaCE.pl genes.gff interproscan.tsv MIN_TPM MIN_EVI_COVERAGE MIN_REF_COVERAGE MAX_AED stringtie_gtf_files > trace.out```
### Arguments:
 1. genes in GFF format
 2. interproscan tab delmited output of Pfam domains
 3. minimum TPM cutoff (0.5)
 4. minimum sample transcript overlap (0.5)
 5. minimum reference transcript overlap (0.5)
 6. maximum annotation edit distance (AED) between sample and reference transcript (0.5)
 7. stringtie .gtf files (one per sample)
 
