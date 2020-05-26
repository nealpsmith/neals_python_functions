import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, Formula, conversion
pandas2ri.activate()
from rpy2.robjects.packages import importr
import numpy as np
deseq = importr('DESeq2')
'''
Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
'''

to_dataframe = robjects.r('function(x) data.frame(x)')

class py_DESeq2 :
    '''
    DESeq2 object through rpy2
    input:
    count_matrix: should be a pandas dataframe with each column as count, and a id column for gene id
        example:
        id    sampleA    sampleB
        geneA    5    1
        geneB    4    5
        geneC    1    2
    design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
                treatment
    sampleA1        A
    sampleA2        A
    sampleB1        B
    sampleB2        B
    design_formula: see DESeq2 manual, example: "~ treatment""
    gene_column: column name of gene id columns, exmplae "id"
    '''
    def __init__(self, count_matrix, design_matrix, design_formula) :

        self.dds = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_count_df = None
        # self.gene_column = self.count_matrix.index
        self.gene_id = count_matrix.index
        self.samplenames = count_matrix.columns
        self.count_matrix = pandas2ri.py2rpy(count_matrix)
        self.design_matrix = pandas2ri.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, 
                                        colData=self.design_matrix,
                                        design=self.design_formula)


    def run_deseq(self, **kwargs):
        self.dds = deseq.DESeq(self.dds, **kwargs)


    def get_deseq_result(self, contrast=None, **kwargs):

        self.comparison = deseq.resultsNames(self.dds)
        if contrast:
            if len(contrast)==3:
                contrast = robjects.numpy2ri.numpy2ri(np.array(contrast)) 
            else:
                assert len(contrast) == 2, 'Contrast must be length of 3 or 2'
                contrast = robjects.ListVector({None:con for con in contrast})
            print('Using contrast: ', contrast)
            self.deseq_result = deseq.results(self.dds, contrast = contrast, **kwargs)
        else:
            self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result = to_dataframe(self.deseq_result)
        self.deseq_result = conversion.rpy2py(self.deseq_result)
        
        return(self.deseq_result)