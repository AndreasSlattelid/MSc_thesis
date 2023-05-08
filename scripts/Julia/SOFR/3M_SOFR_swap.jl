using Revise

#statistics and distributions
using Random 
using Distributions
using Statistics

#data-wrangeling:
using DataFrames

#for numerical integration:
using QuadGK

#plotting histogram and LaTeX labels:
using Plots 
using LaTeXStrings

#---------------------------------------------------------------------
#time parameters
T0 = 1/12
T1 = 4/12
T2 = 7/12
T3 = 10/12
timepoints = [T0, T1, T2, T3]
n_steps = 10000

#We use that Vasicek is ATS, i.e P(t,T) = exp(-A(t,T)-B(t,T)r(t))

function r_Vasicek(alpha, m,sigma, r_t, time_interval, n_steps)
    """
    Args: 
        #Vasicek parameters:
        alpha{Float64}: speed of reversion 
        m{Float64}: long term mean level 
        sigma{Float64}: volatility 
        r_t{Float64}: initial value of r = (r(u))
        
        #time:
        time_interval (vector): the time interval we model measured in years.
        n_steps (int): number of timesteps we partition over
    
    Returns: 
        it simulates the process: r = (r(u)), for u in [t_start, t_end]
    """

    #time partition: 
    t_start = time_interval[1]
    t_end = time_interval[2]
    dt = (t_end-t_start)/n_steps
    n_steps = length(collect(t_start:dt:t_end))

    #initializing r: 
    r = zeros(n_steps)
    r[1] = r_t

    #Standard normal rv's
    Z = rand(Normal(0,1), n_steps)

    for i in 2:n_steps
        r[i] = r[i-1] - alpha*(m-r[i-1])*dt + sigma*sqrt(dt)*Z[i]
    end

    return r
end

function B_ZCB(t,T)
    ans = -(1/alpha)*(exp(-alpha*(T-t))-1)
    return ans
end

function A_ZCB(t,T)
    integral, _ = quadgk(u -> B_ZCB(u,T)^(2), t,T)
    ans = m*B_ZCB(t,T) - m*(T-t) - (1/2)*sigma^(2)*integral 
    return ans
end

function P(t,T, r_t)
    ans = exp(-A_ZCB(t,T) -B_ZCB(t,T)*r_t)
    return ans
end

#--------------------------------------------------------------
#Calculating f^3M(t,S,T), again using ATS structure:
#f^3M(t,S,T) = 1/(T-S)*(exp(A(t,S,T) + B(t,S,T)r(t))-1)
function Sigma1(u,t,S,T)
    ans = exp(-alpha*(S-u)) - exp(-alpha*(T-u))
    return ans   
end

function Sigma2(u,t,S,T)
    ans = 1-exp(-alpha*(T-u))
    return ans
end

function B(t,S,T)
    ans = (1/alpha)*(exp(-alpha*(S-t))-exp(-alpha*(T-t)))
    return ans
end

function A(t,S,T)
    first_part = m*(T-S) - m*B(t,S,T)
    c1_2, _  = quadgk(u -> Sigma1(u, t, t, S)^(2), t,S) 
    c2_2, _  = quadgk(u -> Sigma2(u, t, S, T)^(2), S,T)

    ans = first_part + (1/2)*(sigma^2/alpha^2)*(c1_2 + c2_2)

    return ans
end


function f_3M(t,S,T, r_t)
    ans = (1/(T-S))*(exp(A(t,S,T) + B(t,S,T)*r_t) - 1)
    return ans
end 

#time t-value of kappa in 3M SOFR-futures rate swap
function kappa_t(t, r_t)
    "
    Args: 
        t{Float64}: vector of time points i.e [0,T1]
        r_t{Float64}: vector of realization of interest rate model
    Returns: 
        above = sum(P(t,T_{i}*f^{3M}(t,T_{i-1}, T_{i})), i = 1:n)
        below = sum(P(t,T_{i}), i = 1:3)
        kappa_t_3M_SOFR = above/below
    "
    ZCB_prices = map(T -> P(t,T, r_t), timepoints[2:end])
    f_3M_rates = map((x, y) -> f_3M(t, x, y, r_t), timepoints[1:end-1], timepoints[2:end])
    above = sum(ZCB_prices.*f_3M_rates)
    below = sum(ZCB_prices) 
    ans = above/below
    return ans
end


#Vasicek parameters: 
alpha = 0.25
m = 0.035
r_0 = 0.0425
sigma = 0.02

#time:
t_start = 0 
t_end = T0
dt = (t_end-t_start)/n_steps
t = collect(t_start:dt:t_end)

#initialization of simulation
n_sim = 2

R = zeros(length(t), n_sim) #Vasicek rates
K = zeros(length(t), n_sim) #fixed rate kappa for each pair (t,r_t)

Random.seed!(1234)
for i in 1:n_sim
    #Vasicek realization:
    r = r_Vasicek(alpha, m, sigma, r_0, [t_start,t_end], n_steps)
    #collect the time t rate and time t kappa:
    R[:, i] = r
    K[:, i] = map((x,y)-> kappa_t(x,y), t, r)
end

R
K[1]

#plot of rates
plot(R, layout = (1,1), 
        legend = false,
        title = L"t \mapsto r(t),\alpha = 0.25, m = 0.035, \sigma = 0.02, r_{0} = 0.0425 "
        )
xticks!([0, 10_000/2 ,10_000], ["0", L"\frac{T_{0}}{2}", L"T_{0}"])


#plot of kappa_t
plot(K, layout=(1,1), 
        legend = false, 
        title = L"t \mapsto \kappa_{t}^{3M-SOFR},\alpha = 0.25, m = 0.035, \sigma = 0.02, r_{0} = 0.0425"
        )
xticks!([0, 10_000/2 ,10_000], ["0", L"\frac{T_{0}}{2}", L"T_{0}"]) 


