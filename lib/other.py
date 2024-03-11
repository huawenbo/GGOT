import torch
import pandas as pd
import numpy as np


dataset_list = ['GSE48452', 'GSE2565', 'LUAD', 'XJTUSepsis', 'COAD', 'GSE154918']
for dataset in dataset_list:
    # dataset = 'GSE154918'
    data = pd.read_csv(f'./data/process/{dataset}/{dataset}.csv', index_col=0)
    print(data.shape, end=' ')
    group_index = pd.read_csv(f'./data/source/{dataset}/group.csv', index_col=0)
    result = torch.load(f'./result/{dataset}/result.pt')
    args = result['args']
    def get_group(dataname):
        group_index = pd.read_csv(f'./data/source/{dataname}/group.csv', index_col=0)
        group = []
        group_class = []
        for i in range(len(group_index)):
            group.append(np.arange(group_index.values[i, 0], group_index.values[i, 1]))
            group_class.append(group_index.index[np.ones(group_index.values[i, 1] - group_index.values[i, 0], dtype=int) * i])
        group_class = np.concatenate(group_class).astype(str)
        return group, group_class
    group, group_class = get_group(dataset)
    tipping_point = np.argmax(result['gwd'])
    print(tipping_point, result['lwd'].shape)
    data_sort = data.loc[data.index[np.argsort(np.abs(result['lwd'][tipping_point]))[::-1]]]
    data_sort.to_csv(f'./data_sort/{dataset}_sort.csv')

for dataset in dataset_list:
    a = pd.read_csv(f'./data_sort/{dataset}_sort.csv', index_col=0).index[:200]
    b = pd.read_csv(f'./result/{dataset}/trigger_molecules.csv', index_col=0).index
    print((a == b).all(), end=' ')