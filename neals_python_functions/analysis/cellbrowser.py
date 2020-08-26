import pandas as pd
import os

def _make_cell_browser_files(
	adata,
	output_filepath,
	which_meta = "all",
	cluster_column = "leiden_labels",
	embedding = "umap", 
	which_vars = ["auroc", "mean_logExpr", "mean_logExpr_other", "log_fold_change", "percentage", "percentage_other"],
	de_selection_var = "auroc",
	de_selection_cutoff = 0.5) :
	
	# Make the metadata
	if which_meta == "all" :
		cell_meta = adata.obs
	else :
		cell_meta = adata.obs[which_meta]

	# Make the expression matrix
	expr_mtx = pd.DataFrame(adata.X.toarray(), index = adata.obs_names).T
	expr_mtx["gene"] = adata.var_names

	# Adjust the columns so genes is first
	cols = expr_mtx.columns.tolist()
	cols = cols[-1:] + cols[:-1]
	expr_mtx = expr_mtx[cols]

	# Need to make one of the columns the "cluster" columns
	cell_meta = cell_meta.rename(columns = {cluster_column : "cluster"})

	# make cell_name column and make it first
	cell_meta['cellName'] = cell_meta.index
	cols = cell_meta.columns.tolist()
	cols = cols[-1:] + cols[:-1]
	cell_meta = cell_meta[cols]

	# Now lets make the umap file
	embedding_file = pd.DataFrame(adata.obsm["X_{embedding}".format(embedding = embedding)], columns = ["x", "y"])
	embedding_file["cellName"] = adata.obs_names
	embedding_file = embedding_file[["cellName", "x", "y"]]

	# Now the DE file
	de_df = pd.DataFrame(columns = which_vars + ["gene","cluster"])

	for clust in set(adata.obs[cluster_column]) :
		df_dict = {}

		for var in which_vars :
			df_dict[var] = adata.varm["de_res"]["{var}:{clust}".format(var = var, clust = clust)]

		df = pd.DataFrame(df_dict)
		df["gene"] = adata.var_names
		df["cluster"] = clust
		df = df[df[de_selection_var] > de_selection_cutoff]
		de_df = de_df.append(df, ignore_index = True)
	
	# Adjust so cluster is first column
	# Adjust the columns so genes is first
	de_df = de_df[["cluster", "gene"] + which_vars]

	if not os.path.exists(output_filepath) :
		os.mkdir(output_filepath)
	
	# Write out the files to the path
	expr_mtx.to_csv(os.path.join(output_filepath, "expr_mtx.csv.gz"), compression= "gzip", index = False)
	cell_meta.to_csv(os.path.join(output_filepath, "meta_data.csv"), index=False)
	embedding_file.to_csv(os.path.join(output_filepath, "embedding.csv"), index = False)
	de_df.to_csv(os.path.join(output_filepath, "de_data.csv"), index = False)


def _make_conf(
	output_filepath,
	name,
	priority = 10,
	tags = ["10X"],
	shortLabel = 'cell browser',
	expr_mtx = "expr_mtx.csv.gz",
	geneIdType = "auto",
	meta = "meta_data.csv",
	enumFields = ["cluster"],
	coord_file = "embedding.csv",
	coord_label = "embedding",
	clusterField = "cluster",
	labelField = "cluster",
	marker_file = "de_data.csv",
	marker_label = "",
	showLabels = True,
	radius = 5,
	alpha = 0.3,
	unit = "log2_CPM",
	matrixType = "auto"
	) :
	
	coord_dict = {
	"file" : '{coord_file}'.format(coord_file = coord_file),
	"flipY" : False,
	"shortLabel" : coord_label
	}
	marker_dict = {
	"file" : '{marker_file}'.format(marker_file = marker_file),
	"shortLabel" : marker_label
	}

	with open(os.path.join(output_filepath, "cellbrowser.conf"), "w") as f :
		f.write(
			"""
			name = '{name}'
			priority = {pri}
			tags = {tags}
			shortLabel = '{shortLabel}'
			exprMatrix='{expr_mtx}'
			geneIdType='{gid}'
			meta = '{meta}'
			enumFields = {enum}
			coords = [
				{coord_dict}
			]
			clusterField='{clusterField}'
			labelField='{labelField}'
			markers=[
			{marker_dict}
			]
			showLabels = {showLabels}
			radius = {radius}
			alpha = {alpha}
			unit = '{unit}'
			matrixType = '{matrixType}'

			""".format(name = name, pri = priority, tags = tags, shortLabel = shortLabel, expr_mtx = expr_mtx,
				gid = geneIdType, meta = meta, enum = enumFields, coord_dict = coord_dict,
				clusterField = clusterField, labelField = labelField, marker_dict = marker_dict,
				showLabels=  showLabels, radius = radius, alpha = alpha, unit = unit, matrixType = matrixType)
			)
		f.close()

	# This is dumb, but it works to remove the whitespace that ends up in the text file
	with open(os.path.join(output_filepath, "cellbrowser.conf"), "r") as f :
		lines = f.readlines()
		lines = [line.replace('	', '') for line in lines]
		f.close()

	with open(os.path.join(output_filepath, "cellbrowser.conf"), "w") as f:
		f.writelines("\n".join(lines))
		f.close()

def _make_browser(data_filepath, browser_filepath, run = True) :
	import subprocess
	# subprocess.call(["cd", data_filepath])
	print('Running cbBuild')
	completed_process = subprocess.run(["cbBuild", "-o", browser_filepath], cwd=data_filepath,
									   stdout=subprocess.PIPE,
									   stderr=subprocess.STDOUT)
	print(completed_process.stdout.decode())
	completed_process.check_returncode()

	# Change the files
	_swap_files(browser_filepath)

	if run :
		# Re-run it
		completed_process = subprocess.run(["cbBuild", "-o", browser_filepath, "-p", "1234"], cwd = data_filepath,
										   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		print(completed_process.stdout.decode())
		completed_process.check_returncode()

def _swap_files(browser_filepath) :
	import shutil
	# Lets delete the original files and replace with new ones
	for direc in ["css", "ext", "js"] :
		shutil.rmtree(os.path.join(browser_filepath, direc))

		shutil.copytree(os.path.join(os.path.dirname(__file__), "db", "cell_browser", direc), os.path.join(browser_filepath, direc))


	os.remove(os.path.join(browser_filepath, "index.html"))
	shutil.copy2(os.path.join(os.path.dirname(__file__), "db", "cell_browser", "index.html"), os.path.join(browser_filepath))


# A function that creates Kamil's cellbrowser
def make_kamil_browser(
	adata,
	browser_filepath,
	browser_name = "cell_browser",
	which_meta = "all",
	cluster_column = "leiden_labels",
	embedding = "umap", 
	which_vars = ["auroc", "mean_logExpr", "mean_logExpr_other", "log_fold_change", "percentage", "percentage_other"],
	de_selection_var = "auroc",
	de_selection_cutoff = 0.5,
	run_browser = True,
	**kwargs # To be passed to the _make_conf function
	) :

	### CHECK THE DATA ###
	# Make sure all of the vars are in the anndata
	all_vars = set([name.split(":")[0] for name in adata.varm["de_res"].dtype.names])
	if not all(var in all_vars for var in which_vars) :
		raise ValueError("Not all DEG parameters are in data")

	# Make sure embedding exists
	if not "X_{embedding}".format(embedding = embedding) in adata.obsm :
		raise ValueError("{embedding} is not in the data".format(embedding = embedding))

	# Make sure cluster column exists
	if not cluster_column in adata.obs.columns :
		raise ValueError("{clustcol} is not in obs".format(clustcol = cluster_column))
	
	# Lets make the metadata
	if not which_meta == "all" :
		if not set(which_meta).issubset(adata.obs.columns) :
			raise ValueError("not all metadata in obs")


	# Lets make the cell browser files
	_make_cell_browser_files(adata, output_filepath = browser_filepath, which_meta = which_meta, 
		cluster_column = cluster_column, embedding = embedding, which_vars = which_vars,
		de_selection_var = de_selection_var, de_selection_cutoff = de_selection_cutoff)

	# Lets make the conf file
	_make_conf(output_filepath = browser_filepath, name = browser_name, **kwargs)


	# Create the browser
	_make_browser(data_filepath = browser_filepath, browser_filepath = os.path.join(browser_filepath, "browser"), run = run_browser)