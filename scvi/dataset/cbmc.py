import csv
import os
import numpy as np
import urllib.request
from .dataset import GeneExpressionDataset
import pandas as pd

"""
Without protein information. 
"""

filename = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"


class CbmcDataset(GeneExpressionDataset):
    def __init__(self, type='train'):
        self.save_path = 'data/'
        self.download_name = 'expression.bin'
        self.data_filename = 'expression_%s.npy' % type
        self.labels_filename = 'labels_%s.npy' % type
        self.gene_names = 'genes_names.npy'
        self.download_and_preprocess()
        super(CbmcDataset, self).__init__([np.load(self.save_path + self.data_filename)],
                                            [np.load(self.save_path + self.labels_filename)])

    def download(self):
        url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"
        print("Downloading CBMC data")

        # Create the path to save the data
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        with urllib.request.urlopen(url) as response, open(self.save_path + filename, 'wb') as out_file:
            data = response.read()  # a `bytes` object
            out_file.write(data)

    def preprocess(self):
        print("Preprocessing CBMC data")

        # TODO: Get labels

        expression_data = pd.read_csv(self.save_path + filename, index_col=0, compression='gzip').T
        gene_names = expression_data.columns

        selected = np.std(expression_data.as_matrix(), axis=0).argsort()[-600:][::-1]
        expression_data = expression_data.as_matrix()[:, selected]
        gene_names = gene_names[selected]

        # train test split for log-likelihood scores
        expression_train, expression_test = GeneExpressionDataset.train_test_split(expression_data)

        np.save(self.save_path + "expression_train.npy", expression_train)
        np.save(self.save_path + "expression_test.npy", expression_test)
        # np.save(self.save_path + 'labels_train.npy', c_train)
        # np.save(self.save_path + 'labels_test.npy', c_test)
        np.save(self.save_path + self.gene_names, gene_names)

    def download_and_preprocess(self):
        if not (os.path.exists(self.save_path + self.data_filename) and
                os.path.exists(self.save_path + self.labels_filename)):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()
