using Random 
using Distributions
using Statistics

using LinearAlgebra
using Queryverse
using DataFrames

using Combinatorics
using Plots

function create_array(dims::Array{Tuple{Int,Int},1})
    "
    Args: 
        dims: (array(tuple)), vector of tuples, where each element corresponds to matrix-dimension
    
    Returns:
        Array of matricies A = (M_{1}, …, M_{n})
        Each matrix can take on different dimensions, i.e:
        dimensions = [(m,n), (k,l), (r,q), ...] 
        The array will return matricies of zeros    
    "
    arr = Array{Array{Float64,2}}(undef, length(dims))

    for (i, dim) in enumerate(dims)
        arr[i] = zeros(dim...)
    end

    return(arr)
end

function all_perm(xs, n)
    " 
    Args: 
        xs: (Vector)
        n: (Int), desired length
    
    Returns:
        generates permutation of elements in vector xs of length n
        all_perm([0.0, 1.0], 2) = [0.0, 0.0],[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]
        all_perm([0.0, 1.0], 3) = [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], ...
    "
    return(vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...))))
end 


function OU_CPP(z0::Float64, β::Float64, σ::Float64, λ::Float64, μ::Float64, Δt::Float64,T_end::Float64)
    " 
    Description:
        dZ(t) = -βZ(t)dt + σdW(t) + dI(t)
        I(t): I(t) = ∑_{i=1}^{N(t)}J_{k}, J_{k} ~ Exp(μ), NB Julia parametrize with 1/μ

    Args:    
        z0: (float) inital value of the process Z(t)
        β: (float), mean-retreving parameter
        σ: (float), parameter of Brownian Motion
        λ: (float), jump intensity of process, N(t)~ Pois(λt) 
        Δt: (float) stepsize 
        T_end: (float), for how long the simulation should go
     
    Returns: 
        X(t) = 100exp(-Z(t))
    "

    time = collect(0:Δt:T_end) 
    n = length(time)

    #number of jumps from [0,T_end] on each Δt: N(Δt) ~ Pois(λΔt)
    N = rand(Poisson(λ*Δt),n) 

    #Brownian motion W~N(0,dt) on [0,Δt]
    W = rand(Normal(0,1) ,n)

    #intialising Z(t)
    z = zeros(n)

    z[1] = z0

    #dZ(t) = -βZ(t)dt + σdW(t) + dI(t)
    for i in 2:n 
        ΔI = sum(rand(Exponential(μ), N[i])) - sum(rand(Exponential(μ), N[i-1]))
        z[i] = z[i-1] - β*z[i-1]*Δt + σ*W[i]*Δt + ΔI
    end

    #X(t) = 100exp(-Z(t))
    x = 100*exp.(-z)
    df = DataFrame(time = time, score = x)
    return df
end

function simulation(n_sim, C_ESG, T_end, relevant_times)
    "
    Args:
        n_sim: (int) number of simulations 
        C_ESG: Vector(float) ESG-criteria at time T_{i}            
        relevant_times: (float) vector of relevant times, [T1, T2, T3, ...], percentage of year.
    "

    "
    retruns:
        matrix of wheter or not criteria is meet
        each row in the matrix corresponds to a simulation, i.e
        m = [0,0,0; did not meet any criteria
             0,1,1; met criteria at T2 and T3
             ...  ] 
    "
    
    #store matrix of zeros, row = simulation number, col = agreed observation times
    m = zeros(n_sim, length(relevant_times))
    for i in 1:n_sim
        #df_tmp: general simulation 
        df_tmp = OU_CPP(-log(20/100), -0.05, 0.0,20.0, 1/150.0,1/360, T_end)
        #get the relevant timepoints as df:
        df_relevant = filter(row -> row.time ∈ relevant_times, df_tmp)
        #get the score:
        relevant_score = df_relevant.score
    
        #check if X_{T_{i}} <= C^_{T_{i}}^{ESG} for T_{1}, ..., T_{n}
        ESG_criteria = relevant_score .<= C_ESG 
        ESG_criteria = Float64.(ESG_criteria)
        
        #store ESG_criteria:
        m[i, :] = ESG_criteria
    end

    return(m)
end

function D(i::Int, m::Matrix)
    " 
    Args:
        i: (Int), index in sequence 
        m: (Matrix), matrix containing simulations

    Returns: 
        D(i)-term in in E_{Q}[K_{i}^{ESG}(ω)|F_{t}]    
    "
    if i > size(m)[2]
        return println("You cannot evaluate D outside of agreed contract")
    end
    
    #adjusting for column dimensions in Boolean check:
    m_adj = m[:, 1:i]

    v = all_perm([0.0, 1.0], i)
    possible_patterns = mapreduce(permutedims, vcat, v)
    
    #use the row sum to determine how many errors/fails there are:
    success_sum = collect(0.0:Float64(i))
    
    #p is the indicator of successes for the trial:
    p = zeros(size(possible_patterns)[1], i+1)
    for k in 1:(i+1)
        p[:, k] = Bool[success_sum[k] == sum(possible_patterns[j, :]) for j=1:size(possible_patterns,1)]'
    end
    #turn p into Boolean object so that we can use findall:
    p = Bool.(p)

    #store row dimensions, so that we can initialize array later
    row_dims = zeros(Int, i+1)
    for i in 1:(i+1)
        row_dims[i] = size(possible_patterns[findall(p[:, i]), :])[1]
    end
    
    #initializing the needed dimensions
    dimensions = [(row_dims[k], i) for k in 1:(i+1)] 
    " 
    A: array of matricies, A=(M_{1}, ..., M_{i})
    Let i = 3:
    M_{1}: matrix of patterns giving zero successes  [0,0,0] (1x3)
    M_{2}: matrix of patterns giving one succeses    [0,0,1;
                                                      0,1,0; 
                                                      1,0,0] (3x3)
    M_{3}: matrix of patterns giving two succeses    [1,1,0;
                                                      1,0,1;
                                                      0,1,1] (3x3)                                                
    etc. 
    "
    A = create_array(dimensions) 

    for l in 1:(i+1)
        A[l] = possible_patterns[findall(p[:, l]), :]
    end

    " 
    E_fails: (vector) Expecation of all linear combinations where: 
    E_fails[1]: expectation of all linear combinations giving all fails (1 path)
    E_fails[2]: expectation of all linear combinaition giving fails, but 1 success (multiple paths)
    Let i=3: E_fails[1] = E[fff|F_{t}]
             E_fails[2] = E[ssf|F_{t}] + E[sfs|F_{t}] + E[fss|F_{t}]
    etc.
    "
    E_fails = zeros(i+1)

    for l in 1:(i+1)
        s = 0 
        for j in 1:size(A[l],1)
            s += mean(Bool[A[l][j, :] == m_adj[r, :] for r=1:size(m_adj,1)])
        end
        E_fails[l] = s
    end 
    
    #represents I_{2i}^{Even} = {2,...,2i}
    I_2_Even = collect(2.0:2.0:Float64(2*i))

    #represents vector of ∑_{α∈I_{2i}^{Even}}[i-α], 
    weight = i .- I_2_Even

    #all success, all success but one, all success but two, ... 
    E_success = reverse(E_fails)
    
    #=
    i*E_{Q}[⋂1(A_{l})|F_{t}] + 
    ∑_{α\in I_2_Even}∑_{j_{1} != ... != j_I_α_Even}×
    E_{Q}[(⋂1(A_{l}))^{{j_{1} != ... != j_I_α_Even}}|F_{t}] 
    =#

    ans = i*E_success[1] + sum(weight.*E_success[2:length(E_success)])

    return(ans)
end

relevant_times = [1.25, 2.25, 3.25, 4.25]

C_ESG = [17.6, 16.6, 15.6, 14.55]
m = simulation(10^(6), C_ESG, 4.25, relevant_times)


#----------------------------------------------------------------------------
#κ_{t}^{ESG} and κ: 
#Calculating κ for the above example:

#constants:
d = 0.005 #discount
δ = 1.00  #equidistant distance between T_{i} and T_{i-1}
num_steps = 1:length(relevant_times) #number of relevant steps

#ZCB:
P_T0 = 0.995
P_T1 = 0.985
P_T2 = 0.975
P_T3 = 0.965
P_T4 = 0.955

#vector of ZCB
P = [P_T0, P_T1, P_T2, P_T3, P_T4] 





κ_t_ZCB = (P[1]-P[length(P)])/(δ*sum(P[2:length(P)]))

function κ_t_ESG(i)
    κ_t_ZCB-d*D(i,m)
end

println("(t, κ_t_ZCB, κ_t_ESG, C_ESG, relevant_times)")
for i in 1:length(relevant_times)
    println((i, κ_t_ZCB ,κ_t_ESG(i), C_ESG, relevant_times))
end



#nice parameters:
z0 = -log(20/100)
beta = -0.05
sigma = 0.0
lambda = 20.0
mu = 1/150
dt = 1/360
T_end = 5.0

ESG_index = OU_CPP(z0 , beta, sigma,lambda, mu,dt, T_end)
plot(ESG_index.score)


ESG_plt = ESG_index |> @vlplot(
    title = "X(t)", 
    :line, 
    x = :time, 
    y = {:score, scale = {domain = (12,21)}, 
        field = {"ESG-risk score"}}, 
    width = 375, 
    height = 250 
)

ESG_plt