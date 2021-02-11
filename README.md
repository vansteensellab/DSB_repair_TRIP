# manuscript_DSB_trip

This repository is related to:

*Impact of chromatin context on Cas9-induced DNA double-strand break repair pathway balance* (bioRxiv 2020)

https://www.biorxiv.org/content/10.1101/2020.05.05.078436v1

It has also under revision at Molecular Cell.

**Data availability**

Raw sequencing data: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA686952

Processed data and markdowns: https://osf.io/cywxd/

**Summary**

DNA double-strand break (DSB) repair is mediated by multiple pathways. It is thought that the local chromatin context affects the pathway choice, but the underlying principles are poorly understood. Using a newly developed multiplexed reporter assay in combination with Cas9 cutting, we systematically measured the relative activities of three DSB repair pathways as a function of chromatin context in >1,000 genomic locations. This revealed that non-homologous end-joining (NHEJ) is broadly biased towards euchromatin, while the contribution of microhomology-mediated end-joining (MMEJ) is higher in specific heterochromatin contexts. In H3K27me3-marked heterochromatin, inhibition of the H3K27 methyltransferase EZH2 reverts the balance towards NHEJ. Single-stranded template repair (SSTR), often used for precise CRISPR editing, competes with MMEJ, and is moderately linked to chromatin context. These results provide insight into the impact of chromatin on DSB repair pathway balance, and guidance for the design of Cas9-mediated genome editing experiments.

**Conda environment**


**1. Raw data processing**


**2. Indel & mapping (iPCR) data from the IPRs**

It contains three main scripts that process the data that was created by the DSB_trip_snakamake. 

  1. Parsing_QC
  2. Indel_processing_data
  3. Chromatin_processing_data

The data that was used for these scripts are available in the processed sequencing data folder on https://osf.io/cywxd/
These scripts produce multiple RDS output files that are avaible in the R data output for analysis folder on https://osf.io/cywxd/



**3. Mapping IPRs with Tagmeppr**

**4. Rearrangement detection with tagmentation**

**5. ChIP**

**6. Timeseries**

