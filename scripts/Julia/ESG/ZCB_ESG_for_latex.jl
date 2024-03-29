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
        dims: (array(tuple)), vector of tuples, 
              where each element corresponds to matrix-dimension
    
    Returns:
        Array of matricies A = (M_1, …, M_n)
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


function OU_CPP(z0::Float64, beta::Float64, sigma::Float64, lambda::Float64, mu::Float64, dt::Float64,T_end::Float64)
    " 
    Description:
        dZ(t) = -beta*Z(t)dt + sigma*dW(t) + dI(t)
        I(t): I(t) = sum_{i=1}^{N(t)}J_k, J_k ~ Exp(mu), 
        NB! Julia parametrize with 1/mu

    Args:    
        z0: (float) inital value of the process Z(t)
        beta: (float), mean-retreving parameter
        sigma: (float), parameter of Brownian Motion
        lambda: (float), jump intensity of process, N(t)~ Pois(lambda*t) 
        dt: (float) stepsize 
        T_end: (float), for how long the simulation should go
     
    Returns: 
        X(t) = 100exp(-Z(t))
    "

    time = collect(0:dt:T_end) 
    n = length(time)

    #number of jumps from [0,T_end] on each dt: N(dt) ~ Pois(lambda*dt)
    N = rand(Poisson(lambda*dt),n) 

    #Brownian motion W~N(0,dt) on [0,dt]
    W = rand(Normal(0,1) ,n)

    #intialising Z(t)
    z = zeros(n)

    z[1] = z0

    #dZ(t) = -beta*Z(t)dt + sigma*dW(t) + dI(t)
    for i in 2:n 
        dI = sum(rand(Exponential(mu), N[i])) - sum(rand(Exponential(mu), N[i-1]))
        z[i] = z[i-1] - beta*z[i-1]*dt + sigma*W[i]*dt + dI
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
        relevant_times: (float) vector of relevant times, 
        [T1, T2, T3, ...], percentage of year.
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
        df_relevant = filter(row -> row.time in relevant_times, df_tmp)
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
        D(i)-term in in E_{Q}[K_{i}^{ESG}(omega)|F_{t}]    
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
    A: array of matricies, A=(M_1, ..., M_i)
    Let i = 3:
    M_1: matrix of patterns giving zero successes  [0,0,0] (1x3)
    M_2: matrix of patterns giving one succeses    [0,0,1;
                                                    0,1,0; 
                                                    1,0,0] (3x3)
    M_3: matrix of patterns giving two succeses    [1,1,0;
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
    Let i=3: E_fails[1] = E[fff|F_t]
             E_fails[2] = E[ssf|F_t] + E[sfs|F_t] + E[fss|F_t]
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

    #represents vector of sum_{alpha in I_{2i}^{Even}}[i-alpha], 
    weight = i .- I_2_Even

    #all success, all success but one, all success but two, ... 
    E_success = reverse(E_fails)
    
    #=
    i*E_{Q}[\cap 1(A_{l})|F_t] + 
    sum_[alpha \in I_2_Even]sum_[j_1 != ... != j_I_alpha_Even]x
    E_{Q}[(\cap 1(A_{l}))^[[j_{1} != ... != j_I_alpha_Even]]|F_t] 
    =#

    ans = i*E_success[1] + sum(weight.*E_success[2:length(E_success)])

    return(ans)
end

relevant_times = [1.25, 2.25, 3.25, 4.25]

C_ESG = [17.6, 16.6, 15.6, 14.55]
m = simulation(10^(6), C_ESG, 4.25, relevant_times)


#----------------------------------------------------------------------------
#kappa_{t}^{ESG} and kappa_{t}: 
#Calculating kappa for the above example:

#constants:
d = 0.005 #discount
delta = 1.00  #equidistant distance between T_{i} and T_{i-1}
num_steps = 1:length(relevant_times) #number of relevant steps

#ZCB:
P_T0 = 0.995
P_T1 = 0.985
P_T2 = 0.975
P_T3 = 0.965
P_T4 = 0.955

#vector of ZCB
P = [P_T0, P_T1, P_T2, P_T3, P_T4] 

#ordinary fixed rate kappa form ZCB-swap
kappa_t_ZCB = (P[1]-P[length(P)])/(delta*sum(P[2:length(P)]))

#the ESG-Swap process:
function kappa_t_ESG(i)
    kappa_t_ZCB-d*D(i,m)
end

println("(t, kappa_t_ZCB, kappa_t_ESG, C_ESG, relevant_times)")
for i in 1:length(relevant_times)
    println((i, kappa_t_ZCB ,kappa_t_ESG(i), C_ESG, relevant_times))
end


#meet criteria often:  
C_ESG_2 = [18, 17, 16, 15]
m2 = simulation(10^(5), C_ESG_2, 4.25, relevant_times)

function kappa_t_ESG_2(i)
    kappa_t_ZCB-d*D(i,m2)
end

for i in 1:length(relevant_times)
    println((i, kappa_t_ZCB ,kappa_t_ESG_2(i), C_ESG_2, relevant_times))
end

#does not meet criteria: 
C_ESG_3 = 2.5

m3 = simulation(10^(5), C_ESG_3, 4.25, relevant_times)

function kappa_t_ESG_3(i)
    kappa_t_ZCB-d*D(i,m3)
end

for i in 1:length(relevant_times)
    println((i, kappa_t_ZCB ,kappa_t_ESG_3(i), C_ESG_3, relevant_times))
end

