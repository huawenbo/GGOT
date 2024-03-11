import os.path

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def plot_critical_point(args, wasserstein_distance, x_tricks, size=(6, 2.5), yloc=2.8, path=f'./plots/'):

    # wasserstein_distance = wasserstein_distance * max(wasserstein_distance)
    if not os.path.exists(f'{path}{args.dataset}/'):
        os.makedirs(f'{path}{args.dataset}/')
    plt.figure(figsize=size)
    # plt.tick_params(right=True, top=True)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['font.family'] = 'Arial'
    plt.plot(np.arange(len(wasserstein_distance)), wasserstein_distance, color='r', alpha=0.8)
    plt.scatter(np.arange(len(wasserstein_distance)), wasserstein_distance, color='r', marker='*')

    y_min = min(wasserstein_distance) * 0.6
    y_max = min(wasserstein_distance) * 0.4 + max(wasserstein_distance)
    critical_state = np.argmax(wasserstein_distance)
    pre_state = critical_state - 1
    plt.vlines(critical_state, y_min, wasserstein_distance[critical_state], linestyles='--', colors='green', alpha=0.5)
    plt.vlines(pre_state, y_min, wasserstein_distance[pre_state], linestyles='--', colors='green', alpha=0.5)
    # plt.fill_between([pre_state, critical_state], y_min, y_max, facecolor='red', alpha=0.2, linewidth=0.0)
    plt.fill_between([pre_state, critical_state], y_min,
                     [wasserstein_distance[pre_state], wasserstein_distance[critical_state]], facecolor='red', alpha=0.2, linewidth=0.0)
    # plt.quiver(pre_state + 1.5, y_min * yloc - 0.7, -0.8, 0.2, scale_units='xy', scale=1, width=3e-3, color='green')
    # plt.text(pre_state+2.5, y_min * yloc - 0.8, 'Pre-transition stage', horizontalalignment='right', color='green')
    plt.scatter(critical_state, wasserstein_distance[critical_state], c='b', marker='*', label='Tipping Point', alpha=0.9)
    plt.xticks(np.arange(len(wasserstein_distance)), x_tricks)
    plt.ylim([y_min + 0.3, y_max - 0.3])
    plt.legend()
    plt.xlabel('Stage')
    plt.ylabel('Global Wasser Distance')
    plt.savefig(f'{path}{args.dataset}/gwd.svg', bbox_inches='tight')


if __name__ == '__main__':
    # plot_critical_point(np.array(wass_dist3) / 1e4, ['1', '2', '3'])
    pass
