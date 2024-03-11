import os
import numpy as np
import pandas as pd
from tqdm import tqdm


__all__ = ['GetSubPPI', 'DataProcess']


class GetSubPPI:

    """
    GetSubPPI

    This method is used to get PPI network for Normal Control group
    ----------------------------
    Needs:
     1. ppi_gene_info: The index just select for the high related genes just ppi database, which is restored the
        information.
     2. ppi_matrix: The graph matrix by protein-protein interations networks with a threshold.
     3. ppi_index: The index just select by gene list.
    """

    def __init__(self, dataset, species='Human', gene_list=None, threshold=0.8, expression_ratio=0.9):
        self.dataset = dataset
        self.species = species
        ppi_index_info_source, ppi_index_source = self.get_ppi_index(gene_list, species)
        ppi_matrix = self.get_ppi_matrix(ppi_index_info_source, species)
        index_sub = self.get_sub_ppi(ppi_matrix, threshold, expression_ratio)
        self.ppi_index = ppi_index_source[index_sub]
        self.ppi_gene_info = ppi_index_info_source[index_sub, :]
        self.ppi_matrix = ppi_matrix[index_sub, :][:, index_sub] + np.eye(len(self.ppi_index))
        print('The number of final genes is {}\n'.format(len(self.ppi_index)))

    def __getitem__(self, ind):
        return self.ppi_index[ind]

    def __len__(self):
        return len(self.ppi_index)

    @staticmethod
    def get_ppi_index(gene_list, species):
        print('Starting get ppi index ...')
        map_gene_protein = pd.read_csv(f'./PPI/{species}/map_gene_protein.csv')
        map_gene_protein_full = pd.read_csv(f'./PPI/{species}/map_gene_protein_full.csv')
        gene_full_list = map_gene_protein['gene_id'].values
        index = []
        index_source = []
        for ind, gene_id in enumerate(tqdm(gene_list, ncols=100)):
            if gene_id in gene_full_list:
                index.append([np.where(gene_full_list == gene_id)[0][0], gene_id,
                              map_gene_protein_full.values[np.where(gene_full_list == gene_id)[0][0], 0],
                              map_gene_protein_full.values[np.where(gene_full_list == gene_id)[0][0], 2],
                              map_gene_protein_full.values[np.where(gene_full_list == gene_id)[0][0], 3]])
                index_source.append(ind)
        index_source = np.array(index_source)
        index = np.array(index)
        return index, index_source

    @staticmethod
    def get_ppi_matrix(ppi_index, species):
        print('Starting get ppi matrix ...')
        index = ppi_index[:, 0].astype(int)
        ppi_matrix = np.load(f'./PPI/{species}/ppi_database.npy')
        return ppi_matrix[:, index][index, :]

    def get_sub_ppi(self, ppi_matrix, threshold, expression_ratio):
        index_sub = []
        for ind, vec in enumerate(tqdm(ppi_matrix, ncols=100)):
            vec[ind] = 0
            if np.sum(vec >= threshold) > 0:
                index_sub.append(ind)
        ppi_matrix1 = ppi_matrix[index_sub, :][:, index_sub]
        if not os.path.exists(f'./data/process/{self.dataset}/'):
            os.makedirs(f'./data/process/{self.dataset}/')
        print('Removing the independent nodes with PPI level {0}...'.format(threshold))
        with open(f'./data/process/{self.dataset}/sub_ppi_{expression_ratio}_{threshold}.txt', 'w') as f:
            for ind, vec in enumerate(tqdm(ppi_matrix1, ncols=100)):
                vec[ind] = 0
                f.writelines(str(ind))
                index = np.where(vec >= threshold)[0]
                for ii in index:
                    f.writelines('\t' + str(ii))
                f.writelines('\n')
        return index_sub
    

class DataProcess:

    """
    Provides:
      1. predata: The finally source data without selected ppi, just delete the lacking gene.
      2. ppi_data: The finally data with selected ppi related genes, which is just used by us.
      3. subppi: The class is a dataframe restored the ppi index and ppi networks information and other useful local
         data.
    """

    def __init__(self, dataset, species, is_ppi=False, express_rate=0.9, ppi_level=0.8, express_threshold=0):
        self.express_rate = express_rate
        self.ppi_level = ppi_level
        self.express_threshold = express_threshold
        self.dataset = dataset
        self.species = species

        data = self.data_duplication(dataset=dataset)
        index = self.get_high_expression(data, express_rate, express_threshold)
        data_select, ensemble_id = self.get_predata(data, species, index)

        sub_ppi = GetSubPPI(dataset, species, ensemble_id, ppi_level, express_rate)
        self.data = data_select.iloc[sub_ppi.ppi_index, :]
        self.ppi_graph = sub_ppi.ppi_matrix

        # if is_ppi:
        #     sub_ppi = GetSubPPI(dataset, species, ensemble_id, ppi_level, express_rate)
        #     self.data = data_select.iloc[sub_ppi.ppi_index, :]
        #     self.ppi_graph = sub_ppi.ppi_matrix
        # else:
        #     self.data = data_select
        #     self.ppi_graph = np.ones([len(data_select), len(data_select)])
        self.group, self.stage = self.get_group()
        self.data.to_csv(f'./data/process/{dataset}/{dataset}.csv')
        np.save(f'./data/process/{dataset}/ppi.npy', self.ppi_graph)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.data[index]
    
    @ staticmethod
    def data_duplication(dataset):
        data = pd.read_csv(f'./data/source/{dataset}/{dataset}.csv', index_col=0)
        index_list = pd.unique(data.index)
        if len(data) != len(index_list):
            # print(data.shape, index_list.shape)
            data_mean = data.groupby(data.index).mean()
        else:
            data_mean = data
        return data_mean

    @ staticmethod
    def get_high_expression(data, rate=0.9, express_threshold=0):
        index = []
        n = data.shape[1]
        for i, val in enumerate(data.values):
            length_high_level = len(np.where(val > express_threshold)[0])
            # print(len(np.where(val != 0)[0]))
            if length_high_level / n >= rate:
                index.append(i)
        print('The number of genes with expression level greater than {0} in more than {1} of the {2} samples: {3}'.format(express_threshold, 
                                                                                                                           rate,
                                                                                                                           n,
                                                                                                                           len(index)))
        return index

    @ staticmethod
    def get_predata(data, species, index):
        data_select = data.iloc[index, :]
        # gene_id = pd.read_csv('data/geneids2.csv', index_col=0)  # this is for "select_gene.csv"
        gene_id = pd.read_csv(f'./PPI/{species}/symbol2id.csv', index_col=0)
        gene_select_index = []
        for val in data_select.index:
            if val in gene_id.index:
                gene_select_index.append(val)
        ensemble_id = gene_id['ensembl2'][gene_select_index].tolist()
        data_select = data_select.loc[gene_select_index, :]
        print('The number of genes after mapping the genome is: {0}\n'.format(data_select.shape[0]))
        return data_select, ensemble_id

    def get_group(self):
        group_index = pd.read_csv(f'./data/source/{self.dataset}/group.csv', index_col=0)
        group = []
        for i in range(len(group_index)):
            group.append(np.arange(group_index.values[i, 0], group_index.values[i, 1]))
        stage = group_index.index[1:].values
        return group, stage


if __name__ == '__main__':
    pass
    # sepsis_gene = DataProcess()
