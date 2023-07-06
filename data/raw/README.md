# Counts matrices

## Data origin

In brief: These count matrices come from 12 postmortem humam dorsolateral frontal cortex (Brodmann area 9/46) brain tissue samples. Samples were sequenced using 
one Oxford Nanopore PromethION flow cell per sample, with a median of ~40M mapped reads per sample. The kit used to sequence these samples was PCS111 and flow cell R 9.4.1.
Samples were trimmed and oriented using pychopper. Aligned to the human GRCh38 reference genome using minimap2. Transcript discovery and quantification was performed 
using the Bambu R package. More information about the data and analysis can be found in the methods section of the article associated with this project and 
under the `nextflow_pipeline` directory in the main folder of this GitHub repository.



## File description

`./counts_transcript.txt` - Transcript level counts matrix generated with Bambu using the quantification mode with ENSEMBL annotation 109 (2023)
