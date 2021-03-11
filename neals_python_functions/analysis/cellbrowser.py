import pandas as pd
import numpy as np
import os


def _make_expr_mtx(adata, output_filepath):
    # Make the expression matrix
    expr_mtx = pd.DataFrame.sparse.from_spmatrix(adata.X.T, columns=adata.obs_names)
    expr_mtx["gene"] = adata.var_names

    # Adjust the columns so genes is first
    cols = expr_mtx.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    expr_mtx = expr_mtx[cols]

    expr_mtx.to_csv(os.path.join(output_filepath, "exprMatrix.tsv.gz"), compression="gzip", index=False, sep='\t')


def _make_cell_meta(adata, output_filepath, cluster_column='leiden_labels', which_meta='all'):
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


def _make_embedding(adata, output_filepath, embedding='umap'):
    # Now lets make the embedding file
    embedding_file = pd.DataFrame(adata.obsm[f"X_{embedding}"], columns=["x", "y"])
    embedding_file["cellName"] = adata.obs_names
    embedding_file = embedding_file[["cellName", "x", "y"]]
    embedding_file.to_csv(os.path.join(output_filepath, "embedding.csv"), index=False)


def _make_de_data(
        adata,
        output_filepath,
        which_vars=('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                  'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
        var_info="de_res",
        cluster_column="leiden_labels",
        de_selection_var="auroc",
        de_selection_cutoff=0.5,
        de_selection_top_num=None,
        de_selection_bottom_num=None,
        pval_precision=3,
        round_float=2
):

    # get de columns and filter de_selection_var with de_selection_cutoff
    de_df = pd.DataFrame(columns=list(which_vars) + ["gene", "cluster"])

    for clust in set(adata.obs[cluster_column]):
        df_dict = {}

        for var in which_vars:
            df_dict[var] = adata.varm[var_info]["{clust}:{var}".format(var=var, clust=clust)]

        df = pd.DataFrame(df_dict)
        df["gene"] = adata.var_names
        df["cluster"] = clust
        if de_selection_top_num:
            df = df.sort_values(de_selection_var, ascending=False)[:1000]
        elif de_selection_bottom_num:
            df = df.sort_values(de_selection_var, ascending=True)[:1000]
        else:
            df = df[df[de_selection_var] > de_selection_cutoff]
        de_df = de_df.append(df, ignore_index=True)

    # Adjust so cluster and gene are the first two columns
    de_df = de_df[["cluster", "gene"] + list(which_vars)]

    # make p values display in scientific notation and round other float columns
    pval = 'pseudobulk_p_val'
    for col in which_vars:
        if col != pval:
            de_df[col] = de_df[col].round(round_float)
        elif col == pval:
            de_df[pval] = [np.format_float_scientific(num, precision=pval_precision) for num in de_df[pval]]

    de_df.to_csv(os.path.join(output_filepath, "de_data.csv"), index=False)


def _make_conf(
        output_filepath,
        name,
        priority=10,
        tags="10X",
        shortLabel='cell browser',
        expr_mtx="exprMatrix.tsv.gz",
        geneIdType="auto",
        meta="meta_data.csv",
        enumFields="cluster",
        coord_file="embedding.csv",
        coord_label="embedding",
        clusterField="cluster",
        labelField="cluster",
        marker_file="de_data.csv",
        marker_label="",
        showLabels=True,
        radius=5,
        alpha=0.3,
        unit="log2_CPM",
        matrixType="auto"
):
    coord_dict = {
        "file": '{coord_file}'.format(coord_file=coord_file),
        "flipY": False,
        "shortLabel": coord_label
    }
    marker_dict = {
        "file": '{marker_file}'.format(marker_file=marker_file),
        "shortLabel": marker_label
    }
    tags = [tags]

    with open(os.path.join(output_filepath, "cellbrowser.conf"), "w") as f:
        f.write(
            f"name='{name}'\npriority={priority}\ntags={tags}\nshortLabel='{shortLabel}'\nexprMatrix='{expr_mtx}'\n"
            f"geneIdType='{geneIdType}'\nmeta='{meta}'\nenumFields='{enumFields}'\ncoords=[\n\t{coord_dict}\n]\n"
            f"clusterField='{clusterField}'\nlabelField='{labelField}'\nmarkers=[\n\t{marker_dict}\n]\n"
            f"showLabels={showLabels}\nradius={radius}\nalpha={alpha}\nunit='{unit}'\nmatrixType='{matrixType}'\n"
        )
        f.close()

    # This is dumb, but it works to remove the whitespace that ends up in the text file
    # with open(os.path.join(output_filepath, "cellbrowser.conf"), "r") as f:
    #     lines = f.readlines()
    #     lines = [line.replace('	', '') for line in lines]
    #     f.close()
    #
    # with open(os.path.join(output_filepath, "cellbrowser.conf"), "w") as f:
    #     f.writelines("\n".join(lines))
    #     f.close()


def _make_browser(data_filepath, browser_filepath, run=True):
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


def _swap_files(browser_filepath):
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
def make_kamil_browser(adata,
                       browser_filepath,
                       browser_name="cell_browser",
                       which_meta="all",
                       cluster_column="leiden_labels",
                       embedding="umap",
                       var_info="de_res",
                       which_vars=('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                                   'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
                       de_selection_var="auroc",
                       de_selection_cutoff=0.5,
                       de_selection_top_num=None,
                       de_selection_bottom_num=None,
                       run_browser=True,
                       pval_precision=3,
                       round_float=2,
                       **kwargs):
    prepare_cb_files(adata, browser_filepath, which_meta, cluster_column, embedding, var_info, which_vars,
                     de_selection_var, de_selection_cutoff, de_selection_top_num, de_selection_bottom_num,
                     pval_precision, round_float)
    run_cbBuild(browser_filepath, browser_name, run_browser, **kwargs)


def prepare_cb_files(
        adata,
        browser_filepath,
        which_meta="all",
        cluster_column="leiden_labels",
        embedding="umap",
        var_info="de_res",
        which_vars=('auroc', 'log2Mean', 'log2Mean_other', 'log2FC', 'percentage', 'percentage_other',
                    'percentage_fold_change', 'mwu_U', 'mwu_pval', 'mwu_qval'),
        de_selection_var="auroc",
        de_selection_cutoff=0.5,
        de_selection_top_num=None,
        de_selection_bottom_num=None,
        pval_precision=3,
        round_float=2,
):
    ### CHECK THE DATA ###
    # Make sure all of the vars are in the anndata
    all_vars = set([name.split(":")[1] for name in adata.varm[var_info].dtype.names])
    if not all(var in all_vars for var in which_vars):
        print('If the following output is a list of cluster identifiers, you may need to update Pegasus to 1.0.0 or '
              'higher and rerun de_analysis.')
        print(all_vars)
        raise ValueError(f"Not all {var_info} parameters are in data")

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
    selection_truth_vals = [not de_selection_bottom_num, not de_selection_top_num, not de_selection_cutoff]
    assert sum(selection_truth_vals) <= 1, 'Must pick zero or one of three selection variables to use'

    # create output directory
    os.makedirs(browser_filepath, exist_ok=True)

    # Lets make the cell browser files
    _make_expr_mtx(adata, browser_filepath)
    _make_cell_meta(adata, browser_filepath, cluster_column, which_meta)
    _make_embedding(adata, browser_filepath, embedding)
    _make_de_data(adata, browser_filepath, which_vars, var_info, cluster_column,
                  de_selection_var, de_selection_cutoff, de_selection_top_num, de_selection_bottom_num,
                  pval_precision, round_float)



def run_cbBuild(browser_filepath,
                browser_name="cell_browser",
                run_browser=True,
                **kwargs  # To be passed to the _make_conf function
                ):
    # Lets make the conf file
    _make_conf(output_filepath=browser_filepath, name=browser_name, **kwargs)

    # Create the browser
    _make_browser(data_filepath=browser_filepath, browser_filepath=os.path.join(browser_filepath, "browser"),
                  run=run_browser)
