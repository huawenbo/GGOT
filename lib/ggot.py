import torch
import numpy as np
import os
from multiprocessing.pool import Pool
import pandas as pd


class GGOT:

    def __init__(self, args, data, group, ppi_graph=None):
        self.args = args
        self.data = data
        self.ppi_graph = ppi_graph
        self.group = group
        print('1 Getting the mean and covariance')
        self.mu, self.sigma = self.get_stage_gauss(self.data, group=self.group)
        self.mu = torch.from_numpy(self.mu)
        self.sigma = torch.from_numpy(self.sigma)
        print('2 Calculating the wasserstein distance and optimal transport map')
        self.sigma_0, self.sigma_0_root, _, self.sigma_0_error = self.mat_sqrt(self.sigma[0], 0)
        if not os.path.exists(f'./result/{args.dataset}/'):
            os.makedirs(f'./result/{args.dataset}/')
        state = [[0, i] for i in range(1, len(self.group))]
        if args.multicore:
            pool = Pool(processes=args.n_cores)
            result_multi = pool.map(self.wasserstein_distance, state)
            pool.close()
            pool.join()
        else:
            result_multi = []
            for val in state:
                result_multi.append(self.wasserstein_distance(state))
        self.result = self.result_convert(result_multi)
        self.molecules_trigger = self.get_trigger_molecules(args.dataset, args.trigger_molecules)

    def result_convert(self, result_multi):
        result = {'args': self.args, 'gwd': [], 'lwd': [], 't_map': [], 'error': []}
        for ind, val in enumerate(result_multi):
            result['gwd'].append(val[ind + 1][0].cpu().numpy())
            result['lwd'].append(val[ind + 1][1].cpu().numpy())
            # result['t_map'].append(val[ind + 1][2].cpu().numpy())
            # result['error'].append(val[ind + 1][2].cpu().numpy())
        result['gwd'] = np.array(result['gwd'])
        result['lwd'] = np.array(result['lwd'])
        # torch.save(result, f'./result/{self.args.dataset}/result_distance.pt')
        # result['t_map'] = np.array(result['t_map'])
        # result['error'] = np.array(result['error'])
        torch.save(result, f'./result/{self.args.dataset}/result.pt')
        return result

    def get_gauss(self, data_stage):
        num = data_stage.shape[1]
        mu = data_stage.mean(axis=1).values
        data_stage_loss = data_stage - data_stage.mean(axis=1).values.reshape(-1, 1)

        # sigma = np.zeros([len(reference_data_loss), len(reference_data_loss)])
        # for i in tqdm(range(len(reference_data_loss))):
        #     for j in range(i, len(reference_data_loss)):
        #         # biased estimation
        #         # sigma[i, j]=(reference_data_loss.values[i, :] * reference_data_loss.values[j, :]).mean()
        #         # unbiased estimation
        #         sigma[i, j]=(reference_data_loss.values[i, :] * reference_data_loss.values[j, :]).sum() / (num - 1)
        #         sigma[j, i]=sigma[i, j]

        # sigma = (reference_data_loss.values @ reference_data_loss.values.T) / (num)  # biased estimation
        sigma = (data_stage_loss.values @ data_stage_loss.values.T) / (num-1)  # unbiased estimation
        if self.args.is_ppi:
            sigma = sigma * self.ppi_graph
        return mu, sigma

    @ staticmethod
    def mat_sqrt(mat, index, is_output=True):
        l, v = torch.linalg.eig(mat)
        l_real = torch.real(l)
        v_real = torch.real(v)
        l_real[l_real < 0] = 0
        l_real = l_real + 1e-20
        matrix = v_real @ torch.diag(l_real) @ v_real.T
        matrix_sqrt = v_real @ torch.diag(torch.sqrt(l_real)) @ v_real.T
        matrix_sqrt = (matrix_sqrt + matrix_sqrt.T) / 2
        matrix_inv = v_real @ torch.diag(1 / torch.sqrt(l_real)) @ v_real.T
        matrix_inv = (matrix_inv + matrix_inv.T) / 2
        error = torch.norm(mat - matrix_sqrt @ matrix_sqrt)
        if is_output:
            print('The error of the sqrt of matrix {0} is: {1:.3f}'.format(index, error))
        return matrix, matrix_sqrt, matrix_inv, error
    
    def get_stage_gauss(self, data, group):
        mu = np.zeros([len(group), data.shape[0]])
        sigma = np.zeros([len(group), data.shape[0], data.shape[0]])
        for ind, index_group in enumerate(group):
            mu[ind, :], sigma[ind, :, :] = self.get_gauss(data.iloc[:, index_group])

        if not os.path.exists(f'./gauss/{self.args.dataset}'):
            os.makedirs(f'./gauss/{self.args.dataset}')
        np.save(f'./gauss/{self.args.dataset}/mu.npy', mu)
        np.save(f'./gauss/{self.args.dataset}/sigma.npy', sigma)
        return mu, sigma

    def wasserstein_distance(self, ind):
        ind_0, ind_i = ind[0], ind[1]
        _, sigma_i_source = self.sigma[ind_0], self.sigma[ind_i]
        sigma_i, sigma_i_root, _, _ = self.mat_sqrt(sigma_i_source, index=ind_i)
        _, sigma_case, sigma_case_inv, error_case = self.mat_sqrt(self.sigma_0_root @ sigma_i @ self.sigma_0_root, index=ind_i, is_output=False)
        # sigma_root_case = mat_sqrt(sigma2)
        w2 = torch.trace(self.sigma_0) + torch.trace(sigma_i) - 2 * torch.trace(sigma_case)
        # print(ind_i, torch.trace(self.sigma_0), torch.trace(sigma_i), 2 * torch.trace(sigma_case))
        w2_sub = torch.diag(self.sigma_0) + torch.diag(sigma_i) - 2 * torch.diag(sigma_case)
        return {ind_i: [w2, w2_sub, self.sigma_0_root @ sigma_case_inv @ self.sigma_0_root, error_case]}

    # def wasserstein_distance(self, ind):
    #     ind1, ind2 = ind[0], ind[1]
    #     mu1, sigma1 = self.mu[ind1], self.sigma[ind1]
    #     mu2, sigma2 = self.mu[ind2], self.sigma[ind2]
    #     w1 = torch.norm(mu1-mu2) * torch.norm(mu1-mu2)
    #     sigma_case, sigma_case_inv, error_case = self.mat_sqrt(self.sigma_root @ sigma2 @ self.sigma_root, index=ind2)
    #     # sigma_root_case = mat_sqrt(sigma2)
    #     w2 = torch.trace(sigma1) + torch.trace(sigma2) - 2 * torch.trace(sigma_case)
    #     return {ind2: [w1 + w2, self.sigma_root @ sigma_case_inv @ self.sigma_root, (w1, w2, error_case)]}

    def get_trigger_molecules(self, dataset, num_max):
        tipping_point = np.argmax(self.result['gwd'])
        # lwd_tipping = result['lwd'][tipping_point]
        lwd_tipping = np.abs(self.result['lwd'][tipping_point])
        print(f'The tipping point of dataset {dataset}: {tipping_point + 1}')
        molecules_trigger = pd.DataFrame(self.data.index[np.argsort(lwd_tipping)][-num_max:][::-1].values, columns=['Symbol'])
        molecules_trigger['importance'] = np.arange(1, num_max+1)[::-1]
        molecules_trigger.to_csv(f'./result/{dataset}/trigger_molecules.csv', index=False)
        return molecules_trigger

    def predict_sample(self, group, sample_index):
        pass


class ModelExplanation:

    def __init__(self, args, data, group, ppi_graph=None):
        pass



class CrossPrediction:

    def __init__(self):
        pass


if __name__ == '__main__':
    pass
