import pegasus as pg
import numpy as np
import pandas as pd
import concurrent.futures
from sklearn.metrics.cluster import adjusted_rand_score
import random


# Use Rand index to determine leiden resolution to use
def rand_index_plot(adata, rep = "pca", resamp_perc = 0.9, resolutions = (0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9)) :
    # Copy the anndata object to not write any leiden information to original data
    adata = adata.copy()
    rand_indx_dict = {}
    indx_array = adata.obs.index.values
    n_cells = range(adata.shape[0])
    resamp_size = round(adata.shape[0] * resamp_perc)
    for resolution in resolutions :
        try :
            adata.obs = adata.obs.drop("leiden_labels", axis = 1)
        except ValueError as e:
            pass
        pg.leiden(adata, resolution=resolution, rep=rep)

        with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
            futures = [executor.submit(collect_samples, n_cells, resamp_size, indx_array, adata, resolution, rep) for i in range(25)]
            rand_list = [f.result() for f in futures]

        rand_indx_dict[str(resolution)] = rand_list
        print("Finished {res}".format(res=resolution))
    return rand_indx_dict


def collect_samples(n_cells, resamp_size, indx_array, adata, resolution, rep):
    samp_indx = indx_array[random.sample(n_cells, resamp_size)]
    samp_data = adata[samp_indx]

    true_class = samp_data.obs["leiden_labels"]

    pg.leiden(samp_data, resolution=resolution, rep=rep)
    new_class = samp_data.obs["leiden_labels"]
    return adjusted_rand_score(true_class, new_class)


def myeloid_scores(data) :
    from . import gene_sets
    sets = {"DC1_score" : gene_sets.dc1_genes,
    "DC2_score" : gene_sets.dc2_genes,
    "DC3_score" : gene_sets.dc3_genes,
    "DC4_score" : gene_sets.dc4_genes,
    "DC5_score" : gene_sets.dc5_genes,
    "pDC_score" : gene_sets.pdc_genes,
    "mono1_score" : gene_sets.mono1_genes,
    "mono2_score" : gene_sets.mono2_genes,
    "mono3_score" : gene_sets.mono3_genes,
    "mono4_score" : gene_sets.mono4_genes}

    for key, val in sets.items() :
        data.obs[key] = _score_cells(data, val)

def _score_cells(data, gene_set) :
    # Get rid of genes that aren't in data
    gene_set = [gene for gene in gene_set if gene in data.var_names]
    print(gene_set)
    # Limit the data to just those genes
    dat = data[:,gene_set].X
    dat = dat.toarray()
    mean = dat.mean(axis=0)
    var = dat.var(axis=0)
    std = np.sqrt(var)

    with np.errstate(divide="ignore", invalid="ignore"):
        dat = (dat - mean) / std
    dat[dat < -5] = -5
    dat[dat > 5] = 5

    scores = dat.mean(axis = 1)
    return(scores)