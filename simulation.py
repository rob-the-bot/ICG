# %%
"""
Created on Nov 18, 2024
@author: Robert Wong
"""

import logging
from tempfile import TemporaryFile
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import loadmat, savemat
import pandas as pd
from sklearn.linear_model import LinearRegression

# set up logging
logging.basicConfig(level=logging.INFO)

# %% functions which generate covariance matrices


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
    cov = np.random.normal(loc=loc, scale=0.01, size=(n, n))
    cov = (cov + cov.T) / 2
    # change the diagonals to be 1
    cov[np.diag_indices_from(cov)] = 1
    # assert matrix is positive definite
    assert np.all(np.linalg.eigvals(cov) > 0)

    return cov


# %% generate correlated data
N = 4096

np.random.seed(42)
means = np.random.randn(N)
cov = generate_simple_cov(N, loc=0.02)
X = np.random.multivariate_normal(mean=means, cov=cov, size=5000).T
corr_hat = np.corrcoef(X)
logging.info(f"avg correlation value: {np.triu(corr_hat, k=1).mean():.3f}")

# %% send data to MATLAB for ICG calculation

with TemporaryFile() as f:
    logging.info(f"Writing to {f.name}")
    savemat(f.name, {"X": X.astype(np.float32)})
    with TemporaryFile() as f2:
        logging.info(f"Calculating in MATLAB")
        # call matlab to calculate results
        result = subprocess.run(
            ["matlab", "-batch", f"run_ICG('{f.name}', '{f2.name}')"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logging.error(f"MATLAB error: {result.stderr}")
        else:
            logging.info(f"MATLAB output: {result.stdout}")

        res = loadmat(f2.name)["res"].squeeze()

# %% Calculate summary statistics
df = []
for i, X in enumerate(res):
    df.append((i, np.var(X, axis=1).mean()))
df = pd.DataFrame(df, columns=["level", "averaged std"])
df

# %% fit linear regression and plot

reg = LinearRegression().fit(df["level"].values[..., None], np.log2(df["averaged std"]))

df["log std"] = np.log2(df["averaged std"])
fig, ax = plt.subplots(figsize=(3, 2), layout="constrained")
sns.regplot(data=df, x="level", y="log std", color="black", ax=ax, marker="o")
ax.set(
    xlabel="log2(ensemble size)",
    ylabel="log2(averaged std)",
    title=f"Slope: {reg.coef_[0]:.3f}",
)
sns.despine()
