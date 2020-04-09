import pegasus as pg
import numpy as np
import pandas as pd
# Use Rand index to determine leiden resolution to use
def rand_index_plot(adata, rep = "pca", resamp_perc = 0.9, resolutions  = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]) :
	import random
	from sklearn.metrics.cluster import adjusted_rand_score

	# Copy the anndata object to not write any leiden information to original data
	adata = adata.copy()
	rand_indx_dict = {}
	indx_array = adata.obs.index.values
	n_cells = range(adata.shape[0])
	resamp_size = round(adata.shape[0] * resamp_perc)
	for resolution in resolutions :
		try :
			adata.obs = adata.obs.drop("leiden_labels", axis = 1)
		except :
			pass
		pg.leiden(adata, resolution = resolution, rep = rep)

		rand_list = []
		for iter in range(25) :
			samp_indx = random.sample(n_cells, resamp_size)
			samp_indx = indx_array[samp_indx]
			samp_data = adata[samp_indx]

			true_class = samp_data.obs["leiden_labels"]

			pg.leiden(samp_data, resolution = resolution, rep = rep)
			new_class = samp_data.obs["leiden_labels"]

			rand_list.append(adjusted_rand_score(true_class, new_class))
		rand_indx_dict[str(resolution)] = rand_list
		print("Finished {res}".format(res = resolution))
	return(rand_indx_dict)