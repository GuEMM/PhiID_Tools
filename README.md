# PhiID Tools for discrete time series analysis

Julia implementation of Integrated Information Decomposition framework for discrete multivariate time series.

- Implemented for analysis of two groups of data ordered in columns, requires user to indicate the columns related to first group and second.

- The software is implemented to operate using parallel computation to explore any possible bipartition of each group of data. 

- Required Julia packedges are: Random,Distributions,Combinatorics,CausalityTools, DelimitedFiles.
