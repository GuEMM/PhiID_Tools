# PhiID Tools for discrete time series analysis

Julia implementation of Integrated Information Decomposition framework for discrete multivariate time series.

- Implemented for analysis of two groups of variables ordered in columns in the same data file, requires user to indicate the columns related to first and second group.

- The software operates using parallel computation through Julia threads.

- Required Julia packedges are: Random,Distributions,Combinatorics,CausalityTools, DelimitedFiles.

File IIT_PhiR_tools.jl contain all functions needed to compute the Phi-ID atoms from a multivariate system composed by two different groups of element. For example, Inhibitory and Excitatory groups of neurons as done in Menesse and Torres (2024).

File PHIID_atoms.jl run the analysis for a given data file. User should parse the index of the data file in the data directory.

Analysis return a two files. Containing the following information
    - First row indicating the Time Delayed Mutual Information of the full system (N discrete variables analyzed).
    NxM matrix containing:
    - Columns: 16 [Float32] columns contained the 16 information atoms of PhiID (valid for biparitions), [Int32] size of one of the partitions, 2 [Float32] columns with each partition Entropy.
    - Rows: Data for each bipartition analized.

# EXAMPLE

Presents the julia script for analysis of spiking time series of neurons from an EEG neuronal network model (Menesse y Torres, 2023), see https://github.com/GuEMM/EEG_model.git.

The processing and visualization notebook present scripts for computing the statistics of Phi-ID analysis of each possible bipartition on each neuronal group. Also, include functions to compute measure such as Revised Integrated Information, Information Transfer, Information differentiation, Non-synergistic redundancy and others. See Menesse y Torres, 2023.

The notebook Example_computing_info_measures (in python), shows how to select the information atoms values obtained for each bipartition using the IIT_PhiR_tools.jl, and compute different information measures such as Integrated information theory, storage, non-synergistic redundancy and others.

# REFERENCE

Menesse G, Torres J. 2023. Information dynamics efficiently discriminates high γ-rhythms in EEG brain waves. DOI: https://arxiv.org/abs/2311.13977
