# PhiID Tools for discrete time series analysis

Julia implementation of Integrated Information Decomposition framework for discrete multivariate time series.

- Implemented for analysis of two multi-variate variables ordered in columns, requires user to indicate the columns related to first and second variable.

- The software operates using parallel computation through Julia threads.

- Required Julia packedges are: Random,Distributions,Combinatorics,CausalityTools, DelimitedFiles.

File IIT_PhiR_tools.jl contain all functions needed to compute the Phi-ID atoms from a multivariate system composed by two different groups of element. For example, Inhibitory and Excitatory groups of neurons as done in Menesse and Torres (2024).

File PHIID_atoms.jl run the analysis for a given data file. User should parse the index of the data file in the data directory.

Analysis return a two files. Containing the following information
    - First row indicating the Time Delayed Mutual Information of the full system (N discrete variables analyzed).
    NxM matrix containing:
    - Columns: 16 [Float32] columns contained the 16 information atoms of PhiID (valid for biparitions), [Int32] size of one of the partitions, 2 [Float32] columns with each partition Entropy.
    - Rows: Data for each bipartition analized.
