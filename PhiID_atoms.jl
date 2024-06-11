include("IIT_PhiR_tools.jl");

gridsize = "50x50"

#Directory containing the discrete time series
ddirload = "tau_mu_map_"*gridsize*"/"

allfiles = readdir(ddirload)

Num = parse(Int32,ARGS[1])

file = [i for i in allfiles if occursin("_States.txt", i)][Num]

PHItau = []

Tau = [1,10,100]

TTiE = Int32(200000)
TTiI = Int32(200000)

#Create directory to save PhiID atoms results for each bipartition of excitatory (E) and inhibitory (I) groups
IITdirE = "PhiID_"*gridsize*"_12E_datos_$(TTiE)/"
IITdirI = "PhiID_"*gridsize*"_9I_datos_$(TTiI)/"

isdir(IITdirE) || mkdir(IITdirE)
isdir(IITdirI) || mkdir(IITdirI)

#Matrices for read spiking train time series
XE = zeros(Int32,TTiE,12)
XI = zeros(Int32,TTiI,9)
	
mu = parse(Float32,split(file,"_")[6])

taur = parse(Float32,split(file,"_")[8])

#Number of columns related to time series of the first group (I Group)
NI = 9

#All the remains columns will be analyzed as a second group (E group)

#Measure PhiID using different time delays
for ta in length(Tau):-1:1

	tau = Int32(Tau[ta])
 
        filex = "PhiID_atoms_tau_$(tau)_taur_$(taur)_mu_$(mu).txt"
	
	if isfile(IITdirE*filex) == true

		println("Exists")
	
	else
		println("Compute")

		XE[:] = Array{Int32,2}(readdlm(ddirload*file,',',Int8))[1:TTiE,(NI+1):end]
        
		PHIID_parallel_Atoms!(XE,tau,IITdirE*filex)

		XI[:] = Array{Int32,2}(readdlm(ddirload*file,',',Int8))[1:TTiI,1:NI]

		PHIID_parallel_Atoms!(XI,tau,IITdirI*filex)        
		
	end
end
