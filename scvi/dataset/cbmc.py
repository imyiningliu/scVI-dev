import os
import numpy as np
import urllib.request
import scipy.sparse as sp_sparse
import csv
import gzip

from .dataset import GeneExpressionDataset


class CbmcDataset(GeneExpressionDataset):
    def __init__(self, type='train'):
        self.save_path = 'data/'
        self.download_name = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"
        self.data_filename = 'cbmc_expression_%s.npy' % type
        self.gene_names = 'genes_names.npy'
        self.download_and_preprocess()
        super(CbmcDataset, self).__init__([sp_sparse.csr_matrix(np.load(self.save_path + self.data_filename))],
                                          gene_names=np.load(self.save_path + self.gene_names))

    def download(self):
        url = "https://www.ncbi.nlm.nih.gov/geo/download/" + \
              "?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"
        print("Downloading CBMC data")

        # Create the path to save the data
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        with urllib.request.urlopen(url) as response, open(self.save_path + self.download_name, 'wb') as out_file:
            data = response.read()  # a `bytes` object
            out_file.write(data)

    def preprocess(self):
        print("Preprocessing CBMC data")
        rows = []
        gene_names = []
        with gzip.open(self.save_path + self.download_name, "rt", encoding="utf8") as csvfile:
            data_reader = csv.reader(csvfile, delimiter=',')
            for i, row in enumerate(data_reader):
                rows.append(row[1:])
                gene_names.append(row[0])

        expression_data = np.array(rows[1:], dtype=np.int).T
        gene_names = np.array(gene_names[1:], dtype=np.str)

        selected = np.std(expression_data, axis=0).argsort()[-600:][::-1]
        expression_data = expression_data[:, selected]
        gene_names = gene_names[selected]

        # train test split for log-likelihood scores
        expression_train, expression_test = GeneExpressionDataset.train_test_split(expression_data)

        np.save(self.save_path + "cbmc_expression_train.npy", expression_train)
        np.save(self.save_path + "cbmc_expression_test.npy", expression_test)
        np.save(self.save_path + self.gene_names, gene_names)

    def download_and_preprocess(self):
        if not (os.path.exists(self.save_path + self.data_filename)):
            if not os.path.exists(self.save_path + self.download_name):
                self.download()
            self.preprocess()
