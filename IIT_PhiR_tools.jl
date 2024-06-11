using Random,Distributions,Combinatorics,CausalityTools, DelimitedFiles

function Com!(ID::Array{Int32,1},ss1::Array{Int32,2},ss2::Array{Int32,2},nn::Array{Int32,1},s::Int32,C::Int32)
    
    nn[1] = 0
    
    for ss in combinations(ID, s)
        
        nn[1] += 1
        
        ss1[nn[1],1:s] = ss
        ss2[nn[1],1:C-s] = ID[ID .∉ Ref(ss)]
        
    end
    
    return nothing
end

function varphiR(X,tau,s1,s2)
    
    MI = IMf(X[1:end-tau,:],X[1+tau:end,:])
    
    XX1 = X[:,s1]
    XX2 = X[:,s2]
    
    H1 = ComplexityMeasures.entropy(probabilities(Dataset(XX1)))
    H2 = ComplexityMeasures.entropy(probabilities(Dataset(XX2)))

    Ka = minimum([H1,H2])
    
    MI1 = IMf(XX1[1:end-tau,:],XX1[1+tau:end,:])
    MIc1 = IMf(XX1[1:end-tau,:],XX2[1+tau:end,:])
    MIc2 = IMf(XX2[1:end-tau,:],XX1[1+tau:end,:])
    MI2 = IMf(XX2[1:end-tau,:],XX2[1+tau:end,:])
    
    if Ka>0
        vphi = (MI - MI1 - MI2 + minimum([MIc1,MIc2]))/Ka
    else
        vphi = (MI - MI1 - MI2 + minimum([MIc1,MIc2]))/1e-5 
    end
    
    return vphi
end

function IMf(x,y)
    
    #ps1 = probabilities(Dataset(x));
    
    #h1 = ComplexityMeasures.entropy(ps1)

    #ps2 = probabilities(Dataset(y));

    #h2 = ComplexityMeasures.entropy(ps2)

    #ps3 = probabilities(Dataset([x y]));

    #h3 = ComplexityMeasures.entropy(ps3)

    return ComplexityMeasures.entropy(probabilities(Dataset(x)))+ComplexityMeasures.entropy(probabilities(Dataset(y)))-ComplexityMeasures.entropy(probabilities(Dataset([x y])))
end

function Computs_Caus!(Computs::Array{Float64,1},XX1::Array{Int32,2},XX2::Array{Int32,2},varphiR::Array{Float64,1},MI::Float64,tau::Int32)
    
    Computs[1] = ComplexityMeasures.entropy(probabilities(Dataset(XX1)))
    Computs[2] = ComplexityMeasures.entropy(probabilities(Dataset(XX2)))

    Computs[3] = minimum(Computs[1:2])

    Computs[4] = IMf(XX1[1:end-tau,:],XX1[1+tau:end,:])
    Computs[5] = IMf(XX1[1:end-tau,:],XX2[1+tau:end,:])
    Computs[6] = IMf(XX2[1:end-tau,:],XX1[1+tau:end,:])
    Computs[7] = IMf(XX2[1:end-tau,:],XX2[1+tau:end,:])

    Computs[8] = MI - Computs[4] - Computs[7] + minimum([Computs[5],Computs[6]]) 
    
    if Computs[3]>0
        append!(varphiR,Computs[8]/Computs[3])
    else
        append!(varphiR,Computs[8]/1e-4)
    end
    
    return nothing
end

function MinIJ(X1::Array{Int32},X2::Array{Int32},tau::Int32)
    
    I1 = IMf(X1[1:end-tau,:],X2[tau+1:end,:])
    I2 = IMf(X2[1:end-tau,:],X1[tau+1:end,:])
    
    return minimum([I1,I2])
end

function MinI(X1::Array{Int32},X2::Array{Int32},X::Array{Int32},tau::Int32)
    
    I1 = IMf(X1[1:end-tau,:],X[tau+1:end,:])
    I2 = IMf(X2[1:end-tau,:],X[tau+1:end,:])
    
    return minimum([I1,I2])
end

function MinJ(X1::Array{Int32},X2::Array{Int32},X::Array{Int32},tau::Int32)
    
    I1 = IMf(X[1:end-tau,:],X1[tau+1:end,:])
    I2 = IMf(X[1:end-tau,:],X2[tau+1:end,:])
    
    return minimum([I1,I2])
end

function Compute_Atoms!(X1::Array{Int32},X2::Array{Int32},X::Array{Int32},tau::Int32,B::Array{Float64,1},Atoms::Array{Float64,1})
    
    B[:] = [MinIJ(X1,X2,tau),MinI(X1,X2,X1,tau),MinI(X1,X2,X2,tau),MinI(X1,X2,X,tau),MinJ(X1,X2,X1,tau),
        IMf(X1[1:end-tau,:],X1[1+tau:end,:]),IMf(X1[1:end-tau,:],X2[1+tau:end,:]),IMf(X1[1:end-tau,:],X[1+tau:end,:]),
        MinJ(X1,X2,X2,tau),IMf(X2[1:end-tau,:],X1[1+tau:end,:]),IMf(X2[1:end-tau,:],X2[1+tau:end,:]),IMf(X2[1:end-tau,:],X[1+tau:end,:]),
        MinJ(X1,X2,X,tau),IMf(X[1:end-tau,:],X1[1+tau:end,:]),IMf(X[1:end-tau,:],X2[1+tau:end,:]),IMf(X[1:end-tau,:],X[1+tau:end,:])] 
    
    Atoms[:] = [B[1], -B[1] + B[2], -B[1] + B[3], B[1] - B[2] - B[3] + B[4], -B[1] + B[5], B[1] - B[2] - B[5] + B[6],
        B[1] - B[3] - B[5] + B[7], -B[1] + B[2] + B[3] - B[4] + B[5] - B[6] - B[7] + B[8], -B[1] + B[9], B[1] + B[10] - B[2] - B[9], B[1] + B[11] - B[3] - B[9], -B[1] - B[10] - B[11] + B[12] + B[2] + B[3] - B[4] + B[9], B[1] + B[13] - B[5] - B[9], -B[1] - B[10] - B[13] + B[14] + B[2] + B[5] - B[6] + B[9], -B[1] - B[11] - B[13] + B[15] + B[3] + B[5] - B[7] + B[9], B[1] + B[10] + B[11] - B[12] + B[13] - B[14] - B[15] + B[16] - B[2] - B[3] + B[4] - B[5] + B[6] + B[7] - B[8] - B[9]]    

    return nothing
end

function PHIID_parallel_Atoms!(X::Array{Int32,2},tau::Int32,fileX::String)
    
    allfiles = readdir()

    filesR = [i for i in allfiles if occursin(".txt", i) & occursin("PHIs_", i)]
    
    if length(filesR)>0
        for i in 1:length(filesR)
            rm(filesR[i])
        end
    end
    
    Tt,C = size(X)

    ID = collect(1:C)
    
    # Sizes of bipartition
    Sizes = Array{Int32,2}(reduce(vcat,transpose.(unique([i for i in combinations(append!(ID,ID),2) if (sum(i)==C)&(i[1]-i[2]<=0)]))))
    
    sm = maximum(Sizes[:,1])

    maxS = Int32(binomial(Int32(C),sm))
    
    Ns = Int32(size(Sizes)[1])

    listS = []
    
    Limits = []
    
    for ss in 1:Ns 

        for j in 1:binomial(Int32(C), Int32(ss))

            push!(listS,(ss,j))

        end
    end

    NumProc = size(listS)[1]
    
    Ntr = Threads.nthreads()
    
    if NumProc/Ntr > 1
    	NT = Int(ceil(NumProc/Ntr))
    else
    	while NumProc/Ntr < 1
    		Ntr -=1
    	end
    	NT = Int(ceil(NumProc/Ntr))
    end
    
    for i in 1:Ntr
        if i*NT < NumProc 
            push!(Limits,((i-1)*NT+1,i*NT)) 
        else i*NT >= NumProc
            push!(Limits,((i-1)*NT+1,NumProc)) 
            break
        end
    end

    ID = Array{Int32,1}(collect(1:C))

    Ntrs = size(Limits)[1]
    
    PhiId_Atoms = -1 .*ones(Float64,NumProc,16)
    Ss = zeros(Int32,NumProc)
    K = zeros(Float64,NumProc,2)
    
    Threads.@threads :static for i=1:Ntrs

        idthr = Threads.threadid()

        B = zeros(Float64,16)

        Atoms = zeros(Float64,16);

        sa = zeros(Int32,1)

        ss1 = Array{Int32,2}(zeros(maxS,C))
        ss2 = Array{Int32,2}(zeros(maxS,C))

        for ii in Limits[i][1]:Limits[i][2]

            s = Int32(listS[ii][1])
            j = Int32(listS[ii][2])

            if sa[1]!=s

                nn = zeros(Int32,1)

                Com!(ID,ss1,ss2,nn,s,Int32(C))

                sa .= s

            end

            Compute_Atoms!(X[:,ss1[j,1:s]],X[:,ss2[j,1:(C-s)]],X,tau,B,Atoms)

            PhiId_Atoms[ii,:] = Atoms
            Ss[ii] = s
            K[ii,:] = Float64.([ComplexityMeasures.entropy(probabilities(Dataset(X[:,ss1[j,1:s]]))),ComplexityMeasures.entropy(probabilities(Dataset(X[:,ss2[j,1:(C-s)]])))])

        end
    end
    
    io = open(fileX, "w+")
    
    IM = IMf(X[1:end-tau,1:end],X[1+tau:end,1:end])

    writedlm(io,IM)

    writedlm(io,hcat(PhiId_Atoms,hcat(Ss,K)),",")

    close(io)

    return nothing

end

