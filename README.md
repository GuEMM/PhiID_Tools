# PhiID Tools for discrete time series analysis

Julia implementation of Integrated Information Decomposition framework for discrete multivariate time series.

- Implemented for analysis of two groups of data ordered in columns, requires user to indicate the columns related to first group and second group.

- The software operates using parallel computation through Julia threads.

- Required Julia packedges are: Random,Distributions,Combinatorics,CausalityTools, DelimitedFiles.

File IIT_PhiR_tools.jl contain all functions needed to compute the Phi-ID atoms from a multivariate system composed by two different groups of element. For example, Inhibitory and Excitatory groups of neurons as done in Menesse and Torres (2024).

File PHIID_atoms.jl run the analysis for a given data file. User should parse the index of the data file in the data directory.

Analysis return a two files, each of one with a matrix containing:
  Columns: Time Delayed Mutual Information, 16 columns contained the 16 information atoms of PhiID (only for biparitions), Size of one of the partitions, 2 columns with each partition Entropy.
  Rows: Data for each bipartition analized.
