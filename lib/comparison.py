import os
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.linalg as lg
import copy
import scipy.linalg as slg


def show_network(g, y=np.ones(14) * 0.2, labels=None, pos=None, ax=None, figsize=(5, 5)):
    cmap = plt.cm.coolwarm
    if ax is None:
        plt.figure(figsize=figsize)
        plt.axis('off')
        ax = plt.gca()

    if pos is None:
        pos = nx.random_layout(g)

    nx.draw_networkx_nodes(g, pos, node_size=500, alpha=0.7, node_color=cmap(y), ax=ax, linewidths=0)
    nx.draw_networkx_edges(g, pos, ax=ax, width=0.4)

    if labels is None:
        nx.draw_networkx_labels(g, pos, font_color='black', font_size=15, ax=ax)
    else:
        label_dict = {}
        for i, v in enumerate(g.nodes):
            label_dict[v] = labels[i]
        nx.draw_networkx_labels(g, pos, font_color='black', font_size=15, labels=label_dict, ax=ax)


def wasserstein_distance(mat1, mat2):
    root_1 = slg.sqrtm(mat1)
    # root_2 = slg.sqrtm(mat2)
    dis = np.trace(mat1) + np.trace(mat2) - 2 * np.trace(slg.sqrtm(root_1 @ mat2 @ root_1))
    return dis


class GraphComparison:

    def __init__(self,
                 num_nodes=14,
                 prob=((0.9, 0.1), (0.1, 0.9)),
                 seed=83242):
        self.num_nodes = num_nodes
        self.seed = seed
        self.prob = prob
        self.g0, self.l0, self.g1, self.l1, self.g2, self.l2 = self.get_graph()
        self.plot_graphs()

    def get_graph(self, edge_remove1=((4, 7), (2, 13)), edge_remove2=((1, 3), (7, 11))):
        g0 = nx.stochastic_block_model(sizes=[int(self.num_nodes / 2), self.num_nodes - int(self.num_nodes / 2)],
                                       p=self.prob,
                                       seed=self.seed)
        l0 = nx.laplacian_matrix(g0, range(self.num_nodes)).todense()
        g1 = copy.deepcopy(g0)
        for edge in edge_remove1:
            g1.remove_edge(edge[0], edge[1])
        l1 = nx.laplacian_matrix(g1, range(self.num_nodes))
        l1 = l1.todense()

        g2 = copy.deepcopy(g0)
        for edge in edge_remove2:
            g2.remove_edge(edge[0], edge[1])
        l2 = nx.laplacian_matrix(g2, range(self.num_nodes))
        l2 = l2.todense()
        return g0, l0, g1, l1, g2, l2

    def plot_graphs(self, path='./plots/comparison/'):
        pos = nx.kamada_kawai_layout(self.g0)
        pos[1] = pos[1] - [0.1, 0.1]
        if not os.path.exists(path):
            os.makedirs(path)

        fig1 = plt.figure(figsize=(5, 5))
        plt.axis('off')
        show_network(self.g0, pos=pos, ax=fig1.gca())
        fig1.savefig(f'{path}graph_normal.svg', pad_inches=0, bbox_inches='tight')

        fig2 = plt.figure(figsize=(5, 5))
        plt.axis('off')
        show_network(self.g1, pos=pos, y=np.ones(self.num_nodes) * 0.8, ax=fig2.gca())
        fig2.savefig(f'{path}graph_1.svg', pad_inches=0, bbox_inches='tight')

        fig3 = plt.figure(figsize=(5, 5))
        plt.axis('off')
        show_network(self.g2, pos=pos, y=np.ones(self.num_nodes) * 0.8, ax=fig3.gca())
        fig3.savefig(f'{path}graph_2.svg', pad_inches=0, bbox_inches='tight')

    @staticmethod
    def wasserstein_distance_(mat1, mat2):
        n = len(mat1)
        # adding 1 to zero eigenvalue; does not change results, but is faster and more stable
        l1_tilde = mat1 + np.ones([n, n])
        l2_tilde = mat2 + np.ones([n, n])
        s1_tilde = lg.inv(l1_tilde)
        s2_tilde = lg.inv(l2_tilde)
        root_1 = slg.sqrtm(s1_tilde)
        # root_2 = slg.sqrtm(s2_tilde)
        res = s1_tilde + s2_tilde - 2 * slg.sqrtm(root_1 @ s2_tilde @ root_1)
        return np.trace(res), np.diag(res)


if __name__ == '__main__':
    graph_compare = GraphComparison()
