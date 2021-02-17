import numpy as np
import pandas as pd
import concurrent.futures
from sklearn.metrics.cluster import adjusted_rand_score
import random
import time
import logging
logger = logging.getLogger(__name__)
import leidenalg
import concurrent.futures
import os
from pegasus.tools import construct_graph
from scipy.sparse import csr_matrix


# Use Rand index to determine leiden resolution to use
def rand_index_plot(
        W,  # adata.uns['W_' + rep] or adata.uns['neighbors']
        resamp_perc=0.9,
        resolutions=(0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9),
        max_workers=25,
        n_samples=25,
        random_state=0
    ):
    assert isinstance(W, csr_matrix)
    rand_indx_dict = {}
    n_cells = W.shape[0]
    resamp_size = round(n_cells * resamp_perc)

    for resolution in resolutions:

        true_class = leiden(W, resolution, random_state)

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_collect_samples, W, resolution, n_cells, resamp_size, true_class, random_state)
                       for i in range(n_samples)]
            rand_list = [f.result() for f in futures]

        rand_indx_dict[str(resolution)] = rand_list
        print("Finished {res}".format(res=resolution))
    return rand_indx_dict


def leiden(W, resolution, random_state=0):

    start = time.perf_counter()

    G = construct_graph(W)
    partition_type = leidenalg.RBConfigurationVertexPartition
    partition = leidenalg.find_partition(
        G,
        partition_type,
        seed=random_state,
        weights="weight",
        resolution_parameter=resolution,
        n_iterations=-1,
    )

    labels = np.array([str(x + 1) for x in partition.membership])

    end = time.perf_counter()
    n_clusters = len(set(labels))
    logger.info(f"Finished leiden clustering for res = {resolution}. Get {n_clusters} clusters. "
                f"Time spent = {end-start:.2f}s.")

    return pd.Series(labels)



def _collect_samples(W, resolution, n_cells, resamp_size, true_class, random_state=0):
    samp_indx = random.sample(range(n_cells), resamp_size)
    samp_data = W[samp_indx][:, samp_indx]
    true_class = true_class[samp_indx]
    new_class = leiden(samp_data, resolution, random_state)
    return adjusted_rand_score(true_class, new_class)


def plot_boxplot(dct, figdir, save_name=None):
    import seaborn as sns
    import matplotlib.pyplot as plt

    col1 = []
    col2 = []
    for k in dct.keys():
        for i in dct[k]:
            col1.append(k)
            col2.append(i)

    df = pd.DataFrame({'resolution': col1, 'rand_index': col2})

    sns.boxplot(x='resolution', y='rand_index', data=df)
    plt.axhline(y=0.9, color='r', lw=1.0, linestyle='--')
    plt.xlabel('resolution')
    plt.ylabel('rand_index score')
    plt.savefig(os.path.join(figdir, save_name + '_resolution_boxplot.png'))


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