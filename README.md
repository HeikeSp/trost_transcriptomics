# trost_transcriptomics

## Collection of scripts for the analysis of RNASeq Data

The scripts differ in their input data:

* FPKM values
* expected counts table from rsem per gene


## ``Overview.Rmd``
* uses FPKM values per gene (43034 rows)
	* filtering of genes with few reads: remove rows where maximal value < 1
	* results in 26277 rows remaining
	* do PCA
		* **Observation: PCA of fpkm values without further transformation is not suitable**
		* **Conclusion: Use rld for PCA!**
	* do PCA of log2 values
* uses also expected counts per gene (43034 rows)
* 