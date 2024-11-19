# %%
"""
Created on Nov 18, 2024
@author: Robert Wong
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import loadmat, savemat
import pandas as pd

# %% generate correlated data


def generate_toeplitz(n: int, loc: float = 1) -> np.ndarray:
    # Generate the first row and column for the Toeplitz matrix
    mat = np.nan * np.ones((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = loc ** np.abs(i - j)

    return mat


def generate_simple_cov(n: int, loc: float = 1) -> np.ndarray:
    cov = np.ones((n, n)) * loc
    # change the diagonals to be 1
    cov[np.diag_indices_from(cov)] = 1
    return cov


def generate_random_cov(n: int, loc: float = 1) -> np.ndarray:
    L = np.random.randn(n, n)
    cov = L @ L.T
    return loc * cov


def generate_random_cov2(n: int, loc: float = 1) -> np.ndarray:
    cov = np.random.normal(loc=loc, scale=0.01, size=(n, n))
    cov = (cov + cov.T) / 2
    # change the diagonals to be 1
    cov[np.diag_indices_from(cov)] = 1
    # assert matrix is positive definite
    assert np.all(np.linalg.eigvals(cov) > 0)

    return cov


N = 8196

np.random.seed(42)
means = np.random.randn(N)
# use toeplitz matrix
# cov = generate_toeplitz(N, loc=0.95)
cov = generate_simple_cov(N, loc=0.039)
X = np.random.multivariate_normal(mean=means, cov=cov, size=5000).T
print(f"avg correlation value: {np.triu(np.corrcoef(X), k=1).mean():.3f}")
savemat("input.mat", {"X": X.astype(np.float32)})

# %% Load data

data = loadmat("matlab.mat")["ans"].squeeze()

# Calculate summary statistics and save
res = []
for i, X in enumerate(data):
    res.append((i, np.var(X, axis=1).mean()))

# %% linear regression

df = pd.DataFrame(res, columns=["level", "averaged std"])
df

# %% fit linear regression to find the slope
from sklearn.linear_model import LinearRegression

x = df["level"].values.reshape(-1, 1)
y = np.log2(df["averaged std"].values)
reg = LinearRegression().fit(x, y)
print(f"R^2: {reg.score(x, y)}")
print(f"Slope: {reg.coef_[0]}")

# %% plots

df["log std"] = np.log2(df["averaged std"])
fig, ax = plt.subplots(figsize=(3, 2), layout="constrained")
sns.regplot(data=df, x="level", y="log std", color="black", ax=ax, marker="o")
# set y-scale to log
# ax.set_yscale("log")
sns.despine()

# %%
