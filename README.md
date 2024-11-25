# Groves-CellSys2022
Code and figures related to Groves et al. Cell Systems 2022

## Repository Structure:
- `data`
    - `bulk-rna-seq` : contains all data related to bulk RNA seq analyses, used and generated in `notebooks/bulk` and `notebooks/ParTI-code`.
    - `single-cell-rna-seq`
- `mazebox` : this can be used like a python package and contains several functions used in  `notebooks/single-cell`.
    - `archetypes`
    - `plasticity`
    - `plotting`
    - `preprocessing`
    - `__init__.py`
    - `_settings.py`
- `notebooks`
    - `out` : output from `notebooks`, including a Batch Correction report for combining bulk RNA seq from cell lines and tumors. 
    - `ParTI-code` : MATLAB code to run archetype analysis on bulk RNA seq data. 
        - `human-cell-lines`  
            - `out` : contains original five archetypes shown in paper, as well as other numbers of archetypes tested.
            - `Human-cell-line-RNA-seq.mlx` : Code to generate archetypes for human cell line bulk RNA-seq. 
            - `params.txt` : table showing parameters for each ParTI run, including barcode used in file structure in `out`
        - `thomas-tumors` : contains code and results from ParTI for different numbers of archetypes on bulk RNA-seq of 81 human tumors. Same file structure as `human-cell-lines`.
        - `combined-data` : contains code and results from ParTI for different numbers of archetypes on combined bulk RNA-seq dataset from `notebooks/bulk/Cell-line-tumor-batch-correction-and-clustering.Rmd`, using the file `data/bulk-rna-seq/parti-input/CCLE_Minna_Thomas_COMBAT_no_names.csv`. Same file structure as `human-cell-lines`. 
    - `bulk`
        - `int` : contains intermediary files that were saved for antiquity but are not used in any scripts.
        - `Bulk-Archetype-Tasks.ipynb` : plots related to the inferred tasks for each archetype. This analysis shows **N is enriched in axonogenesis genes and migration, A may correspond to TPCs, and Y corresponds to transit-amplifying PNECs.**
        - `SCLC RNA-seq batch correction-CCLE-Minna.Rmd` : batch corrects CCLE and Minna data, as used in the manuscript.
        - `Cell-line-tumor-batch-correction-and-clustering.Rmd` : batch corrects CCLE, Minna, and Tumor datasets. This analysis is used for ParTI on cell lines and tumors combined.
        - `Compare_cell_lines_tumors_archetypes.Rmd` : comparing cell line archetypes to combined dataset archetypes.
        - `Thomas-Tumors-Bulk-Archetypes.ipynb` : comparing the bulk RNA-seq from 81 human tumors to the original five archetypes in the manuscript (from cell line data). This analysis shows that **tumor samples are encompassed within the cell line archetype space, and tumor data alone shows more variation in ribosomal proteins than cell lines which may not be phenotypically relevant.** 
    - `single-cell`
        - `Human_cell_lines.ipynb`
        - `RPM_tumors.ipynb`
        - `RPM_tumors-CR.ipynb`: using CellRank on RPM time series and tumor data
- `ParTI` : ParTI package code, which is used in `notebooks/ParTI-code`.
- `R_code` : R code from Darren Tyson
- `environment.yml` : conda environment 


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/QuLab-VU/Groves-CellSys2021/sarah)
