Code supporting the manuscript: 

*Multiscale organisation of neuronal activity unifies scale-dependent theories of brain function*

**How to Iteratively coarse-grain your data:**

Run the code "ICG.m" on your time series (rows - variables, columns - time). 

ICG will iteratively coarse-grain the time series by pairing variables by their Pearson correlation, then summing and repeating the procedure on the coarse-grained variables. 

It will output the coarse-grained variables (*activityICG*) across as many levels as possible and the original index of the variables (*outPairID*).

**Extra code**
There is also additional code in the github for running the analysis (e.g., info theory, kurtosis, variance, timescales), network models, and various null models) described in the paper.

**Questions**
If you have any questions, please reach out! 

brandon (dot) munn @sydney.edu.au
