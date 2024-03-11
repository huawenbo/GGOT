import os
import torch
import scipy
import time
import pandas as pd
import argparse
from lib.ggot import *
from lib.dataset import *
from lib.plotting import *
from tqdm import tqdm
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser('GGOT for DCP')

parser.add_argument('-d', '--dataset', type=str, default='GSE154918',
                    help='Dataset to load,  Available: GSE48452, GSE2565, LUAD, COAD, XJTUSepsis, GSE154918')
parser.add_argument('-s', '--species', type=str, default='Human',
                    help='The available database for PPI network: Human, Mus')
parser.add_argument('--is_preprocess', action='store_false')
parser.add_argument('--is_ppi', action='store_false')
parser.add_argument('--express_rate', type=float, default=0.9)
parser.add_argument('--ppi_level', type=float, default=0.8)
parser.add_argument('--express_threshold', type=float, default=0)
parser.add_argument('--trigger_molecules', type=int, default=200)


parser.add_argument('--multicore', action='store_false')
parser.add_argument('--n_cores', type=int, default=10)
parser.add_argument('--cuda', action='store_true')
parser.add_argument('--cross_validation', action='store_true')
# parser.add_argument('--')

args = parser.parse_args()

if __name__ == '__main__':

    time_start = time.time()

    print(f'The GGOT analysis for dataset: {args.dataset}, species: {args.species}\n')
    if args.is_preprocess:
        data_process = DataProcess(dataset=args.dataset,
                                   species=args.species,
                                   is_ppi=args.is_ppi,
                                   express_rate=args.express_rate,
                                   ppi_level=args.ppi_level,
                                   express_threshold=args.express_threshold)

    print('Starting to construct the GGOT model')
    path_result = f'./result/{args.dataset}/result.pt'

    # if not os.path.isfile(path_result):
    #     data_ot = GGOT(args=args,
    #                    data=data_process.data,
    #                    group=data_process.group,
    #                    ppi_graph=data_process.ppi_graph)
    #     result = data_ot.result
    # else:
    #     result = torch.load(f'./result/{args.dataset}/result.pt')
    #     print(result['gwd'])

    data_ot = GGOT(args=args,
                      data=data_process.data,
                      group=data_process.group,
                      ppi_graph=data_process.ppi_graph)
    result = data_ot.result

    print('Starting to analysis the result')
    plot_critical_point(args, result['gwd'], x_tricks=data_process.stage)
    time_end = time.time()
    print(f'This experiment takes {int((time_end - time_start) / 60)} minutes')
