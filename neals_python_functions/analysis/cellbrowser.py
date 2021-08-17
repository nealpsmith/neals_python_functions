import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import warnings
import scanpy as sc
from typing import Iterable, Union


def _make_expr_mtx(
        adata: sc.AnnData,
        output_filepath: str,
):
    filename = os.path.join(output_filepath, "exprMatrix.tsv")

    try:
        import counts_to_csv as ctc
        ctc.counts_to_csv(adata, delimiter='tab', column_orient='obs-names', outfile=filename)
        os.system(f'pigz -f -v {filename}')
    except ModuleNotFoundError as e:
        print('Module counts_to_csv not found. Install for faster TSV creation!')
        print('https://github.com/swemeshy/counts_to_csv')
        # Kamil's write_tsv function from villani-lab/covid/make-cellbrowser.py
        if adata.X.getformat() == 'csr':
            X = adata.X.tocsc()
        else:
            X = adata.X
        f = open(filename, 'w')
        head = ['gene'] + adata.obs.index.values.tolist()
        f.write('\t'.join(head))
        f.write('\n')
        for i in tqdm(range(X.shape[1])):
            f.write(adata.var.index.values[i])
            f.write('\t')
            row = X[:, i].todense()
            row.tofile(f, sep="\t", format="%.7g")
            f.write('\n')
        f.close()
        cmd = f'pigz -f -v {filename}'
        os.system(cmd)


def _make_cell_meta(
        adata: sc.AnnData,
        output_filepath,
        cluster_column: str = 'leiden_labels',
        which_meta: Union[str, Iterable[str]] = 'all'
):
    # Make the metadata
    if which_meta == "all":
        cell_meta = adata.obs
    else:
        cell_meta = adata.obs[which_meta]

    # Need to make one of the columns the "cluster" columns
    cell_meta = cell_meta.rename(columns={cluster_column: "cluster"})

    # make cell_name column and make it first
    cell_meta['cellName'] = cell_meta.index
    cols = cell_meta.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    cell_meta = cell_meta[cols]

    cell_meta.to_csv(os.path.join(output_filepath, "meta_data.csv"), index=False)


def _make_embedding(
        adata: sc.AnnData,
        output_filepath: str,
        embedding: str = 'umap'
):
    # Now lets make the embedding file
    embedding_file = pd.DataFrame(adata.obsm[f"X_{embedding}"], columns=["x", "y"])
    embedding_file["cellName"] = adata.obs_names
    embedding_file = embedding_file[["cellName", "x", "y"]]
    embedding_file.to_csv(os.path.join(output_filepath, "embedding.csv"), index=False)


def _make_de_data(
        adata: sc.AnnData,
        output_filepath: str,
        which_vars: Iterable[str] = ('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                  'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
        var_info: str ="de_res",
        cluster_column: str = "leiden_labels",
        de_selection_var: str = "auroc",
        de_selection_cutoff: Union[int, float] = 0.5,
        de_selection_top_num: int = None,
        de_selection_bottom_num: int = None,
        pval_precision: int = 3,
        round_float: int = 2
):

    # get de columns and filter de_selection_var with de_selection_cutoff
    de_df = pd.DataFrame(columns=list(which_vars) + ["gene", "cluster"])

    for clust in set(adata.obs[cluster_column]):
        df_dict = {}

        try:
            for var in which_vars:
                df_dict[var] = adata.varm[var_info]["{clust}:{var}".format(var=var, clust=clust)]
        except ValueError:
            for var in which_vars:
                df_dict[var] = adata.varm[var_info]["{var}:{clust}".format(var=var, clust=clust)]

        df = pd.DataFrame(df_dict)
        df["gene"] = adata.var_names
        df["cluster"] = clust
        if de_selection_top_num:
            df = df.sort_values(de_selection_var, ascending=False)[:de_selection_top_num]
        elif de_selection_bottom_num:
            df = df.sort_values(de_selection_var, ascending=True)[:de_selection_bottom_num]
        else:
            df = df[df[de_selection_var] > de_selection_cutoff]
        de_df = de_df.append(df, ignore_index=True)

    # Adjust so cluster and gene are the first two columns
    de_df = de_df[["cluster", "gene"] + list(which_vars)]

    # make p values display in scientific notation and round other float columns
    pval_labels = ['P_value', 'pval', 'pvalue', 'P.Value', 'adj.P.Val', 'p_val', 'pVal', 'Chisq_P', 'fdr', 'FDR']
    for col in which_vars:
        if any([p in col for p in pval_labels]):
            de_df[col] = de_df[col].round(round_float)
        else:
            de_df[col] = [np.format_float_scientific(num, precision=pval_precision) for num in de_df[col]]

    de_df.to_csv(os.path.join(output_filepath, "de_data.csv"), index=False)


def _make_conf(
        output_filepath: str,
        name: str,
        priority: int = 10,
        tags: str = "10X",
        shortLabel: str = 'cell browser',
        expr_mtx: str = "exprMatrix.tsv.gz",
        geneIdType: str = "auto",
        meta: str = "meta_data.csv",
        enumFields: str = "cluster",
        coord_file: str = "embedding.csv",
        coord_label: str = "embedding",
        clusterField: str = "cluster",
        labelField: str = "cluster",
        marker_file: str = "de_data.csv",
        marker_label: str = "",
        showLabels: bool = True,
        radius: int = 5,
        alpha: float = 0.3,
        unit: str = "log2_CPM",
        matrixType: str = "auto"
):
    coord_dict = [{
        "file": '{coord_file}'.format(coord_file=coord_file),
        "flipY": False,
        "shortLabel": coord_label
    }]
    marker_dict = [
        {"file": f'{marker_file}', "shortLabel": marker_label},
    ]
    tags = [tags]

    with open(os.path.join(output_filepath, "cellbrowser.conf"), "w") as f:
        f.write(
            f"name='{name}'\npriority={priority}\ntags={tags}\nshortLabel='{shortLabel}'\nexprMatrix='{expr_mtx}'\n"
            f"geneIdType='{geneIdType}'\nmeta='{meta}'\nenumFields='{enumFields}'\ncoords=\t{coord_dict}\n"
            f"clusterField='{clusterField}'\nlabelField='{labelField}'\nmarkers=\t{marker_dict}\n"
            f"showLabels={showLabels}\nradius={radius}\nalpha={alpha}\nunit='{unit}'\nmatrixType='{matrixType}'\n"
        )
        f.close()


def _make_browser(
        data_filepath: str,
        browser_filepath: str,
        run: bool = True
):
    import subprocess
    import sys

    # Need to add path environment variable for subprocess to work
    path = os.environ["PATH"]

    # make sure conda path is included
    conda_path = "/".join([sys.prefix, "bin"])
    path = ":".join([path, conda_path])

    print('Running cbBuild')
    completed_process = subprocess.run(["cbBuild", "-o", browser_filepath], cwd=data_filepath,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT,
                                       env={"PATH": path})
    print(completed_process.stdout.decode())
    completed_process.check_returncode()

    # Change the files
    _swap_files(browser_filepath)

    if run:
        # Re-run it
        completed_process = subprocess.run(["cbBuild", "-o", browser_filepath, "-p", "1234"], cwd=data_filepath,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT,
                                           env={"PATH": path})
        print(completed_process.stdout.decode())
        completed_process.check_returncode()


def _swap_files(browser_filepath: str):
    import shutil
    # Lets delete the original files and replace with new ones
    for direc in ["css", "ext", "js"]:
        shutil.rmtree(os.path.join(browser_filepath, direc))

        shutil.copytree(os.path.join(os.path.dirname(__file__), "db", "cell_browser", direc),
                        os.path.join(browser_filepath, direc))

    os.remove(os.path.join(browser_filepath, "index.html"))
    shutil.copy2(os.path.join(os.path.dirname(__file__), "db", "cell_browser", "index.html"),
                 os.path.join(browser_filepath))


# A function that creates Kamil's cellbrowser
def make_kamil_browser(adata: sc.AnnData,
                       browser_filepath: str,
                       browser_name: str = "cell_browser",
                       which_meta: Union[str, Iterable[str]] = "all",
                       cluster_column: str = "leiden_labels",
                       embedding: str = "umap",
                       var_info: str = "de_res",
                       which_vars: Iterable[str] = ('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                                   'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
                       de_selection_var: str = "auroc",
                       de_selection_cutoff: Union[int, float] = 0.5,
                       de_selection_top_num: int = None,
                       de_selection_bottom_num: int = None,
                       run_browser: bool = False,
                       pval_precision: int = 3,
                       round_float: int = 2,
                       **kwargs):
    prepare_cb_files(adata, browser_filepath, which_meta, cluster_column, embedding,
                     var_info, which_vars, de_selection_var, de_selection_cutoff, de_selection_top_num,
                     de_selection_bottom_num, pval_precision, round_float)
    run_cbBuild(browser_filepath, browser_name, run_browser, **kwargs)


def check_args(
        adata: sc.AnnData,
        browser_filepath: str,
        which_meta: str or iter,
        cluster_column: str,
        embedding: str,
        var_info: str,
        which_vars: Iterable[str],
        de_selection_var: str,
        de_selection_cutoff: Union[int, float],
        de_selection_top_num: int,
        de_selection_bottom_num: int,
        pval_precision: int,
        round_float: int,
):
    ### CHECK THE DATA ###
    # Make sure all of the vars are in the anndata
    all_vars = set([name.split(":")[1] for name in adata.varm[var_info].dtype.names])
    if all([i in adata.obs[cluster_column].cat.categories.values for i in all_vars]):
        warnings.warn('The following output is a list of cluster identifiers. You should to update Pegasus '
                      'to 1.0.0 or higher and rerun de_analysis. For now, it is fine.')
        print(all_vars)
        all_vars = set([name.split(":")[0] for name in adata.varm[var_info].dtype.names])
    if not all(var in all_vars for var in which_vars):
        raise ValueError(f"Not all {var_info} parameters are in data: {all_vars}")

    assert de_selection_var in all_vars
    assert os.path.exists(browser_filepath)

    # Make sure embedding exists
    if not "X_{embedding}".format(embedding=embedding) in adata.obsm:
        raise ValueError("{embedding} is not in the data".format(embedding=embedding))

    # Make sure cluster column exists
    if not cluster_column in adata.obs.columns:
        raise ValueError("{clustcol} is not in obs".format(clustcol=cluster_column))

    # Lets make the metadata
    if not which_meta == "all":
        if not set(which_meta).issubset(adata.obs.columns):
            raise ValueError("not all metadata in obs")

    # make sure de_selection top, bottom, cutoff are correctly set
    selection_truth_vals = [de_selection_bottom_num, de_selection_top_num, de_selection_cutoff]
    assert sum(map(bool, selection_truth_vals)) <= 1, 'Must pick zero or one of three selection variables to use'


def prepare_cb_files(
        adata: sc.AnnData,
        browser_filepath: str,
        which_meta: str = "all",
        cluster_column: str = "leiden_labels",
        embedding: str = "umap",
        var_info: str = "de_res",
        which_vars: Iterable[str] = ('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                    'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
        de_selection_var: str = "auroc",
        de_selection_cutoff: Union[int, float] = 0.5,
        de_selection_top_num: int = None,
        de_selection_bottom_num: int = None,
        pval_precision: int = 3,
        round_float: int = 2,
):
    check_args(adata, browser_filepath, which_meta, cluster_column, embedding, var_info,
               which_vars, de_selection_var, de_selection_cutoff, de_selection_top_num, de_selection_bottom_num,
               pval_precision, round_float)

    os.makedirs(browser_filepath, exist_ok=True)

    # Lets make the cell browser files
    _make_expr_mtx(adata, browser_filepath)
    _make_cell_meta(adata, browser_filepath, cluster_column, which_meta)
    _make_embedding(adata, browser_filepath, embedding)
    _make_de_data(adata, browser_filepath, which_vars, var_info, cluster_column,
                  de_selection_var, de_selection_cutoff, de_selection_top_num, de_selection_bottom_num,
                  pval_precision, round_float)


def run_cbBuild(browser_filepath: str,
                browser_name: str = "cell_browser",
                run_browser: bool = False,
                **kwargs  # To be passed to the _make_conf function
                ):
    # Lets make the conf file
    _make_conf(output_filepath=browser_filepath, name=browser_name, **kwargs)

    # Create the browser
    _make_browser(data_filepath=browser_filepath, browser_filepath=os.path.join(browser_filepath, "browser"),
                  run=run_browser)
