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
        self.download_and_preprocess()
        super(CbmcDataset, self).__init__([np.load(self.save_path + self.data_filename)],
                                            [np.load(self.save_path + self.labels_filename)])

    def download(self):
        # Generating samples according to a ZINB process
        url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"
        # Create the path to save the data
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        with urllib.request.urlopen(url) as response, open(self.save_path + filename, 'wb') as out_file:
            data = response.read()  # a `bytes` object
            out_file.write(data)

# TODO: Unzip the file
#        f = gzip.open(self.save_path + expression, 'rb')
#        file_content = f.read()
#        f.close()

    def preprocess(self):
        print("Preprocessing CBMC data")
        expression = pd.read_csv(self.save_path + filename, index_col=0).T
        gene_names = expression.columns
        selected = np.std(expression.as_matrix(), axis=0).argsort()[-600:][::-1]
        expression = expression.as_matrix()[:, selected]


        expression_train, expression_test, c_train, c_test, n_train, n_test = train_test_split(expression, \
                                                                                               protein.as_matrix(), \
                                                                                               norm_protein.as_matrix(),
                                                                                               random_state=0)

        # train test split for log-likelihood scores
        expression_train, expression_test, c_train, c_test = GeneExpressionDataset.train_test_split(expression_data,
                                                                                                    labels)

        np.savetxt(self.save_path + "expression_train", expression_train)
        np.savetxt(self.save_path + "expression_test", expression_test)
        np.savetxt(self.save_path + "n_train", n_train)
        np.savetxt(self.save_path + "n_test", n_test)
        np.savetxt(self.save_path + "c_train", c_train)
        np.savetxt(self.save_path + "c_test", c_test)

    def download_and_preprocess(self):
        if not (os.path.exists(self.save_path + self.data_filename) and
                os.path.exists(self.save_path + self.labels_filename)):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()
