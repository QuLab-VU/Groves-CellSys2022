# This code is taken from all over the place. It has been regenerated many times, but the final version that is copied here is taken from
# single-cell-cline-plasticity_full_102120.ipynb, which is in mazebox_legacy > mazebox > examples > notebooks > bulk.
import sys
import pandas as pd
import os.path as op

# make sure to use position 1
sys.path.insert(1, "/Users/smgroves/Documents/GitHub/Groves-CellSys2021/")
import mazebox as mb

genes_df = pd.read_csv(
    "../ParTI-code/human-cell-lines/out/5/geneEnrichment_continuous_significant.csv",
    header=0,
)
# genes_df = pd.read_csv('/Users/smgroves/Documents/MATLAB/particode/parti_smg/out/Bulk_final_genes__continuous_significant.csv', header = 0)
genes_df = genes_df.loc[
    genes_df["Significant after Benjamini-Hochberg correction?"] == 1
]
# May also want to subset genes_df by "Is first bin maximal?" == 1

# This depends on the archetypes found by ParTI, and which subtypes they each correspond to
subtype_dict = {1: "SCLC-Y", 2: "SCLC-P", 3: "SCLC-N", 4: "SCLC-A2", 5: "SCLC-A"}

# The function generate_arc_sig_df replaces code found in single-cell-cline-plasticity_full_102120.ipynb
arc_sig_df = mb.ar.generate_arc_sig_df(
    genes_df,
    features_name="Feature Name",
    logfc_name="Mean Difference",
    p_adj_name="P value (Mann-Whitney)",
    subtype_name="archetype #",
    subtype_dict=subtype_dict,
)

geneNames = pd.read_csv(
    "../../data/bulk-rna-seq/parti-input/geneNames_wo_lowgenes.csv",
    header=None,
    index_col=0,
)
archetypes = pd.read_csv(
    "../ParTI-code/human-cell-lines/out/5/322234/322234_arc_full.csv", header=None
)
archetypes = archetypes.T
archetypes.index = geneNames.index

archetypes.columns = ["SCLC-Y", "SCLC-P", "SCLC-N", "SCLC-A2", "SCLC-A"]


sig_matrix_pure, cond, n_genes = mb.ar.generate_signature(
    archetypes,
    arc_sig_df,
    norm=None,
    range_min=20,
    range_max=100,
    range_step=1,
    sort_by_ascending=False,
    print_c=False,
)


sig_matrix_pure.to_csv("../out/bulk/sig_matrix_ParTI_2021.csv")
