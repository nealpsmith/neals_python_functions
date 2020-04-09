import scanpy as sc

def basic_marker_plots(adata, save_folder = ".", save_name = "") :
	sc.settings.figdir = save_folder
	# Make some feature plots
	basic_markers= ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','NKG7','KLRB1','KLRC1','KLRD1','CCL5','PRF1','IFNG','GNLY','TIGIT',
	'LAG3','ENTPD1','HMGB2','HMGN2','STMN1','MKI67','CD69','COTL1','CD4','FOXP3','CTLA4','CD79A','CD79B','CD19','SELL','MS4A1',
	'XBP1','MZB1','IGJ','IGLL5','IGLL1','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-DQB1','HLA-DQA1','HLA-DPB1','HLA-DMA','HLA-DMB','CD22',
	'LYZ','CLEC10A','CD68','C1QA','C1QB','C1QC','KIT','CD69','TPSAB1','CPA3','PIGR','VIM','COL1A2','COL3A1','SPARC','SPARCL1','ENG', 'GZMM',
	'ITGA4', 'CD3G', 'ITGAE', 'ITGAL', 'FMNL1', 'FLNA', 'SELM', 'EPCAM']

	basic_markers = [gene for gene in basic_markers if gene in adata.var_names]

	Cytokine_Interleukin = ['B2M','BMP2','BMP7','CD40LG','FASLG','FLT3LG','IDO1','IDO2','IFNA1','IFNA2','IFNA5','IFNA6','IFNA8','IFNB1','IFNE',
	'IFNG','IFNK','IFNW1','IL10','IL11','IL12A','IL12B','IL13','IL15','IL16','IL17A','IL17B','IL17C','IL17D','IL17F','IL18','IL18BP','IL19','IL1A',
	'IL1B','IL1F10','IL1RN','IL2','IL20','IL21','IL22','IL23A','IL24','IL26','IL27','IL28A','IL29','IL33','IL36A','IL36B','IL36G','IL36RN','IL37',
	'IL4','IL4I1','IL6','IL7','IL8','IL9','ILF3','KITLG','LIF','LTB','MIF','OSM','TGFA','TGFB1','TGFB2','TGFB3','TGFBI','TNF','TNFAIP6','TNFSF10',
	'TNFSF11','TNFSF12','TNFSF13B','TNFSF14','TNFSF15','TNFSF18','TNFSF4','TNFSF8','TNFSF9','CD70','IFNA10','IFNA13','IFNA14','IFNA16','IFNA17',
	'IFNA21','IFNA4','IFNA7','IL25','IL28B','IL3','IL31','IL32','IL34','IL5','ILF2','KLRB1','LTA','TNFSF13']
	Cytokine_Interleukin = [gene for gene in Cytokine_Interleukin if gene in adata.var_names]

	Cytokine_receptor = ['EDA2R','EDAR','EGFR','EPOR','FAS','FLT1','FLT3','FLT4','GHR','IFNAR1','IFNAR2','IFNGR1','IFNGR2','IL10RA','IL10RB','IL11RA',
	'IL11RB','IL12RB1','IL12RB2','IL13RA1','IL13RA2','IL15RA','IL17RA''IL17RB','IL17RC','IL17RD','IL17RE','IL17REL','IL18R1','IL18RAP','IL1R1','IL1R2',
	'IL1RAP','IL1RL1','IL1RN','IL20RA','IL20RB','IL21R','IL22RA1','IL22RA2','IL23R','IL27RA','IL28RA','IL2RA','IL2RB','IL2RG','IL31RA','IL3RA','IL4R',
	'IL5RA','IL6R','IL6ST','IL7R','IL9R','KDR','KIT','LEPR','LIFR','LTBR','MET','MPL','NGFR','OSMR','PDGFRA','PDGFRB','PLEKHO2','PRLR','RELT','TGFBR1',
	'TGFBR2','TGFBR3','TNFRSF10A','TNFRSF10B','TNFRSF10C','TNFRSF10D','TNFRSF11A','TNFRSF11B','TNFRSF12A','TNFRSF13B','TNFRSF13C','TNFRSF14','TNFRSF17',
	'TNFRSF18','TNFRSF19','TNFRSF1A','TNFRSF1B','TNFRSF21','TNFRSF25','TNFRSF4','TNFRSF6B','TNFRSF8','TNFRSF9','ACVR1','ACVR1B','ACVR2A','ACVR2B',
	'AMHR2','BMPR1A','BMPR1B','BMPR2','CCBP2']
	Cytokine_receptor = [gene for gene in Cytokine_receptor if gene in adata.var_names]

	Chemokine = ['CCL1','CCL11','CCL13','CCL16','CCL17','CCL18','CCL19','CCL2','CCL20','CCL21','CCL22','CCL23','CCL24','CCL25','CCL26','CCL27','CCL28','CCL4',
	'CCL5','CCL7','CCL8','CX3CL1','CXCL1','CXCL10','CXCL11','CXCL12','CXCL13','CXCL14','CXCL16','CXCL17','CXCL2','CXCL3','CXCL5','CXCL6','CXCL9','XCL1','PF4',
	'PF4V1','PPBP','IL8','CCL23','XCL2','CCL3','CCL14','CCL15','CCL3L1','CCL3L3','CCL4L1','CCL4L2','XCL2','XCR1','CCR1','CCR10','CCR2','CCR3','CCR4','CCR5',
	'CCR6','CCR7','CCR8','CCR9','CCRL1','CCRL2','CD27','CD40','CD44','CD74','CMKLR1','CNTFR','CRLF2','CSF1R','CSF2RA','CSF2RB','CSF3R','CX3CR1','CXCR1',
	'CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','CXCR7']
	Chemokine = [gene for gene in Chemokine if gene in adata.var_names]

	Inflammasome = ['CIITA','NAIP','NOD1','NOD2','NLRC3','NLRC4','NLRC5','NLRP1','NLRP2','NLRP3','NLRP4','NLRP5','NLRP6','NLRP7','NLRP8','NLRP9','NLRP10',
	'NLRP11','NLRP12','NLRP13','NLRP14','NLRX1','IPAF','AIM2','CASP1','PYCARD','NLRC4','CASP5','NLRP3','NAIP','NLRP1','NLRX1','MEFV','PSTPIP1','PYDC1','CARD8']
	Inflammasome = [gene for gene in Inflammasome if gene in adata.var_names]

	sc.pl.umap(adata, color = basic_markers, save = "_{save_name}_basic_markers.png".format(save_name = save_name))
	sc.pl.umap(adata, color = Cytokine_Interleukin, save = "_{save_name}_cytokine_interleukin.png".format(save_name = save_name))
	sc.pl.umap(adata, color = Cytokine_receptor, save = "_{save_name}_cytokine_receptor.png".format(save_name = save_name))
	sc.pl.umap(adata, color = Chemokine, save = "_{save_name}_chemokine.png".format(save_name = save_name))
	sc.pl.umap(adata, color = Inflammasome, save = "_{save_name}_inflammasome.png".format(save_name = save_name))
