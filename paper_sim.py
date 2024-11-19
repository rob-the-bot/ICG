"""
Created on Nov 15, 2024
@author: Tianfeng Lu
"""

import numpy as np
import matplotlib.pyplot as plt

# X1 = np.random.poisson(size=(2048, 10000)) 
# X2 = np.random.poisson(lam=2, size=(2048, 10000)) # N * T
N = 4096
X = np.random.exponential(size=(N, 5000))# N * T

# %% use mean activity
def compute_var(X):
    corr_mat = np.corrcoef(X)
    corr_mat2 = np.full(corr_mat.shape, -np.inf)
    corr_mat2[np.tril_indices_from(corr_mat, k=-1)] = corr_mat[np.tril_indices_from(corr_mat, k=-1)]
    # hierachy clustering based on the corr_mat
    # coarse-graining
    # add the summed activity
    n_neuron = corr_mat.shape[0]
    assembly_neuron_ids = []
    while n_neuron > 0:
        # find the most correlated pair
        neuron1, neuron2 = np.unravel_index(np.argmax(corr_mat2), corr_mat2.shape)
        # remove the two neurons from the neuron_list
        # add the two neurons to the assembly_neuron_ids
        n_neuron -= 2
        assembly_neuron_ids.append([neuron1, neuron2])
        # update the corr_mat
        corr_mat2[neuron1, :] = -np.inf
        corr_mat2[neuron2, :] = -np.inf
        corr_mat2[:, neuron1] = -np.inf
        corr_mat2[:, neuron2] = -np.inf
    return assembly_neuron_ids

def plot_var(assembly_neuron_ids, X):
    plt.figure()
    assembly_activity = []
    for i in assembly_neuron_ids:
        assembly_activity.append(X[i, :].mean(axis=0))
    assembly_activity = np.array(assembly_activity)
    corr_mat = np.corrcoef(assembly_activity)
    arr = corr_mat[np.tril_indices_from(corr_mat, k=-1)]
    plt.hist(arr, bins=np.arange(-1, 1, 0.01))[-1]
    plt.xlim(-1, 1)
    return assembly_activity, arr.mean(), arr.var()
# %%
assembly_activity = X
n_hierarchies = 11
mean_activity_corr_mean, mean_activity_corr_var = np.zeros(n_hierarchies), np.zeros(n_hierarchies)
for i in range(n_hierarchies):
    assembly_neuron_ids = compute_var(assembly_activity)
    assembly_activity, mean_activity_corr_mean[i], mean_activity_corr_var[i] = plot_var(assembly_neuron_ids, assembly_activity)

# %%
plt.figure(dpi=256)
plt.plot(np.arange(n_hierarchies), mean_activity_corr_mean/mean_activity_corr_mean[0], marker='o')
plt.yscale('log')
plt.xlabel('hierarchies')
plt.ylabel('log(mean_corr/mean_corr[0])')
# %% use sum assembly activity
def plot_var2(assembly_neuron_ids, X):
    plt.figure()
    assembly_activity = []
    for i in assembly_neuron_ids:
        assembly_activity.append(X[i, :].sum(axis=0))
    assembly_activity = np.array(assembly_activity)
    corr_mat = np.corrcoef(assembly_activity)
    arr = corr_mat[np.tril_indices_from(corr_mat, k=-1)]
    plt.hist(arr, bins=np.arange(-1, 1, 0.01))[-1]
    plt.xlim(-1, 1)
    return assembly_activity, arr.mean(), arr.var()
# %%
assembly_activity = X
n_hierarchies = 11
sum_activity_corr_mean, sum_activity_corr_var = np.zeros(n_hierarchies), np.zeros(n_hierarchies)
for i in range(n_hierarchies):
    assembly_neuron_ids = compute_var(assembly_activity)
    assembly_activity, sum_activity_corr_mean[i], sum_activity_corr_var[i] = plot_var2(assembly_neuron_ids, assembly_activity)

# %%
plt.figure(dpi=256)
plt.plot(np.arange(n_hierarchies), sum_activity_corr_var/sum_activity_corr_var[0], marker='o')
plt.yscale('log')
plt.xlabel('hierarchies')
plt.ylabel('log(var_corr/var_corr[0])')