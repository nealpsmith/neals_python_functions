import pandas as pd 
import os
def anndata_from_mtx_files(matrix_csv, meta_csv) :
	'''This function takes in 2 csvs.  The first one corresponds to the location of the matrices.  
	This should be 2 columns: sample: sample names, path:filepath.  The metadata csv at minimum needs one column: sample: sample name.
	All additional columns should correspond to metadata that will end up added to the .obs of the anndata object'''
	import scipy.io as sp 
	import scipy.sparse as sparse
	import anndata

	# Read in the sample info, metadata
	sample_info = pd.read_csv(matrix_csv, index_col = 0)
	meta_data = pd.read_csv(meta_csv, index_col = 0)

	# Make sure all samples in sample info are in metadata
	print(sample_info.index)
	print(meta_data.index)

	try :
		all(meta_data.index == sample_info.index)
	except :
		raise ValueError("matrix csv and metadata csv do not have same samples in them!")
		
	mtx_list = []
	barcode_dict = {}
	feature_dict = {}

	for samp in sample_info.index :
		path = sample_info.loc[samp]["filepath"]
		barcodes = pd.read_csv(os.path.join(path, "barcodes.tsv.gz"),
							   sep="\n", index_col=0, header=None)
		# Add id column
		barcodes["sample"] = [samp for i in range(barcodes.shape[0])]
		barcodes.index = ["_".join([samp, code]) for code in barcodes.index]

		# Add in other metadata
		meta = meta_data.loc[samp]
		for info in meta.index :
			barcodes[info] = [meta[info] for i in range(barcodes.shape[0])]

		# Get the gene name info
		features = pd.read_csv(os.path.join(path, "features.tsv.gz"),
							   sep="\t",
							   header=None, index_col=1, names=["gene_ids", "gene_names", "type"])
		features = features.drop(columns="type")
		
		# Get the actual data matrix
		mtx = sp.mmread(os.path.join(path, "matrix.mtx.gz"))

		# Add everything to the dicts
		barcode_dict[samp] = barcodes
		feature_dict[samp] = features
		mtx_list.append(mtx)

	try :
	    assert len(set(df.shape for df in feature_dict.values())) == 1, "Not all of these samples have the same number of genes in expression matrix"
	except ASsertionError as msg :
	    print(msg)

	# Get the matrices together, transpose for proper shape
	final_mtx = sparse.hstack(mtx_list).T
	final_obs = pd.concat(barcode_dict.values())
	final_var = feature_dict[samp]

	adata = anndata.AnnData(X = final_mtx, obs=  final_obs, var = final_var)

	# Make the genes unique
	adata.var_names_make_unique()

	# print(adata.var)
	return(adata)

