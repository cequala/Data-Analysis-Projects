# Gastric Cancer Samples' with reported CNV, and the well-known major chromosome number for each CNV:
(The postfix 'B' for the 'Sample ID' means 'Baseline', which means the timepoint free from any treatment.)
Sample ID | reported CNV amp. | major chromosome #
A09B | RIKTOR | 5
A10B | CCND3, VEGFA | 6
A14B | CCNE1 | 19
A17B | MDM2 | 12
A18B | FGFR2 | 10
A20B | SMAD3 | 15
A26B | MET | 7
A28B | MDM2 | 12

# File Description 
- infercnv_Baseline_CNV.png
InferCNV(https://github.com/broadinstitute/infercnv) outcome for all samples above. Macrophage from all samples was chosen as the reference.
In the 'Observations (Cells)' part below, expressions over 1 in chr6 for `malilgnent_A10B`, chr12 for (`malilgnent_A17B`, `malilgnent_A28B`) are present as red.
I.e., With the given data, only ('CCND3', 'VEGFA') from the `malilgnent_A10B`, 'MDM2' from the (`malilgnent_A17B`, `malilgnent_A28B`) amplifications are proved via the InferCNV.

- A26B_Numbat_bulk_clones_final.png
Numbat(https://github.com/kharchenkolab/numbat) outcome for the sample 'A26B'.
Unlike the result from the InferCNV, amplification of the CNV in chromosome 7 is observed.

In conclusion, the TCGA subtypes of most samples seem to be GS, rather than CIN.