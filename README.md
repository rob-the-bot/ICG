## Code supporting the manuscript

*Multiscale organisation of neuronal activity unifies scale-dependent theories of brain function*

## How to Iteratively coarse-grain your data

Run the code [`ICG.m`](ICG.m) on your time series (rows - variables, columns - time).
ICG will iteratively coarse-grain the time series by pairing variables by their Pearson correlation,
then summing and repeating the procedure on the coarse-grained variables.
It will output the coarse-grained variables (*activityICG*) across as many levels as possible and the original index of the variables (*outPairID*).

## Simulations

In [`simulation.py`](simulation.py) and calling MATLAB externally, we show the scaling exponent 1.5 can be recreated by random activity with weakly correlated covariance matrix.
