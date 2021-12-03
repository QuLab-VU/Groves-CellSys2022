import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os.path as op
import numpy as np
from lifelines.statistics import logrank_test

survival_type = "progression-free_survival"
tumor_dir = "../data/bulk-rna-seq/"
tumor_df = pd.read_csv(op.join(tumor_dir, "thomas_tumors.csv"), header=0, index_col=0)
survival = pd.read_csv(
    op.join(tumor_dir, "response_clusters.csv"), header=0, index_col=0
)
survival = survival.reindex(tumor_df.columns)
print(survival.head())

#%%

survival = survival.dropna(subset=[survival_type])
tumor_df = tumor_df[survival.index]
print(tumor_df.head())

from sklearn.mixture import GaussianMixture


def binarize(data):
    gm = GaussianMixture(n_components=2)
    for gene in data.columns:
        d = data[gene].values.reshape(data.shape[0], 1)
        gm.fit(d)

        # Figure out which cluster is ON
        idx = 0
        if gm.means_[0][0] < gm.means_[1][0]:
            idx = 1

        data[gene] = gm.predict_proba(d)[:, idx]
    return data


from lifelines import KaplanMeierFitter

f = np.vectorize(lambda x: 0 if x < threshold else 1)
from lifelines.plotting import add_at_risk_counts
from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter

## Example Data
# durations = [5,6,6,2.5,4,4] <- replace with overall survival
# event_observed = [1, 0, 0, 1, 1, 1] <- replace with binary high or low expression
out_dir = "./figures/"
for g in ["IFRD1", "MYC", "YAP1"]:
    # plt.hist(tumor_df.loc[g])
    # plt.show()
    durations = survival[survival_type].values
    # binarize expression
    # expr = binarize(tumor_df.T)
    threshold = tumor_df.loc[g].median()
    expr = f(tumor_df.loc[g])
    df = pd.DataFrame(expr, index=survival.index, columns=["e"])
    ## create a kmf object

    ## Fit the data into the model
    groups = df["e"]
    i1 = groups == 0  ## group i1 , having the pandas series  for the 1st cohort
    i2 = groups == 1  ## group i2 , having the pandas series  for the 2nd cohort
    plt.figure()
    ## fit the model for 1st cohort
    kmf_low = KaplanMeierFitter()

    kmf_low.fit(durations[i1], [1] * len(durations[i1]), label=f"Low {g}")
    a1 = kmf_low.plot(color="b")
    ## fit the model for 2nd cohort
    kmf_high = KaplanMeierFitter()

    kmf_high.fit(durations[i2], [1] * len(durations[i2]), label=f"High {g}")
    kmf_high.plot(color="r", ax=a1)
    plt.xlabel("Time (Months)", fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    plt.xticks(rotation=0, fontsize=12)
    plt.ylim((0, 1))
    add_at_risk_counts(kmf_low, kmf_high, ax=a1)
    plt.hlines(y=0.5, xmin=plt.xlim()[0], xmax=plt.xlim()[1], linestyle="--")
    plt.vlines(x=np.median(durations[i1]), ymin=plt.ylim()[0], ymax=0.5, linestyle="--")
    plt.vlines(x=np.median(durations[i2]), ymin=plt.ylim()[0], ymax=0.5, linestyle="--")

    r = logrank_test(
        durations[i1],
        durations[i2],
        event_observed_A=[1] * len(durations[i1]),
        event_observed_B=[1] * len(durations[i2]),
    )
    plt.yticks(rotation=0, fontsize=10)
    plt.xticks(rotation=0, fontsize=10)
    plt.xlabel("At Risk", fontsize=10)

    plt.text(
        plt.xlim()[0],
        1.01 * plt.ylim()[1],
        f"p-value = {np.round(r.p_value, 4)}",
        fontsize=12,
    )
    plt.tight_layout()
    if survival_type == "progression-free_survival":
        plt.title("Progression-free Survival for Patients")
        plt.savefig(op.join(out_dir, f"{g}_pf_survival_curve.pdf"))
    elif survival_type == "overall_survival":
        plt.title("Overall Survival for Patients")
        plt.savefig(op.join(out_dir, f"{g}_survival_curve.pdf"))
