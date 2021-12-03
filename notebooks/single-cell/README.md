# Single Cell RNA Sequencing Analyses

## Human Cell Lines

I split the human cell line analysis into multiple notebooks. Each notebook saves an `h5ad` version of `adata` to be opened in the next notebook.

Key points for ordering analyses:
1. Archetype analysis applied directly to single cell data (and PCA transformation to bulk PCA space):
    1. We want to use MAGIC imputation (which can affect the X layer) 
    2. Possibly scanorama batch correction (which definitely affects X) 
    3. Possibly cell cycle regression (also definitely affects X)
2. Archetype analysis based on the bulk signatures (signature subtyping):
    1. The data needs to be preprocessed in the same way as the bulk data. This means filtered, normalized, and log-transformed, which is done with Dropkick. 
    2. We do not want to use MAGIC imputed or batch corrected values for scoring. 
    2. If we run RNA velocity prior to this, mazebox will transform both the data and an extrapolated timepoint into archetype (signature) space. 
    4. We may want to consider using some sort of filtering of genes to see if the single cell data is in the same space as the bulk data.
3. Velocity analysis and CellRank
    1. The velocity analysis is using spliced and unspliced counts, but it also takes into account neighborhood graphs, which depend on X. 
    2. Therefore, we want to apply this to each sample independently (by batch), and apply it to the neighborhood calculation before MAGIC or batch correction.
    3. CellRank runs moments, recover_dynamics from scvelo, then dynamical velocity. Recover_dynamics uses velocity genes estimated from the steady state model (implicitly, or you can run steady state velocity first) *or* a list of genes, which could be all of them.
    4. CellRank then requires additional kernels, such as pseudotime, to be calculated, which may have their own restrictions on preprocessing.
4. Plasticity analysis
    1. Should be run after CellRank, but maybe should only use velocity calculations (like how I did it previously). Alternatively, I may be able to incorporate more information into generating a transition matrix, which could be utilized to strengthen the plasticity calculations.