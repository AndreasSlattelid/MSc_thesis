using Random 
using Distributions
using Statistics
using StatsPlots #qqplot

#for matrix operations and linear programming
using LinearAlgebra
using JuMP  #lp-problem setup
using HiGHS #lp-solver

#data-wrangeling:
using DataFrames

#for numerical integration:
using QuadGK

#rerun calculations easier:
using Revise

#plotting histogram and LaTeX labels:
using Plots 
using LaTeXStrings

#--------------------------------------------------------------------------#
#Vasicek parameters: 
alpha = 0.25
m = 0.035
r_t = 0.0425
sigma = 0.02

#time parameters
t = 0
S = 1/12
T1M = 2/12
T2M = 3/12
T = 4/12

function Sigma1(u,t,S,T)
    ans = exp(-alpha*(S-u)) - exp(-alpha*(T-u))
    return ans   
end

function Sigma2(u,t,S,T)
    ans = 1-exp(-alpha*(T-u))
    return ans
end

function int_r_start_stop(low,up, t)
    "
    The integral: integral(r(u)du, low, up) as described in Eq (5.5) p.66
    Args: 
        low{Float64}: lower integration limit
        up{Float64}:  upper integtation limit
    Returns: 
        the integral: int_low_up r(u)du
    "
    if t > low
        return "Please chose t <=low"
    end

    #Sigma1 is N(0, int_t_low c1_2 du), Sigma2 is N(0, int_low_up c2_2 du)
    c1_2, _  = quadgk(u -> Sigma1(u, t, low,up)^(2), t,low) 
    c2_2, _  = quadgk(u -> Sigma2(u, t, low,up)^(2), low,up)
    
    c1 = sqrt(c1_2)
    c2 = sqrt(c2_2)
    Z = rand(Normal(0,1))
    ans = ((r_t-m)/alpha)*(exp(-alpha*(low-t))-exp(-alpha*(up-t))) + m*(up-low) + sigma/alpha*(c1*Z + c2*Z)
    return ans
end

function integrand_E_Q_r(u, r_t, t)
    ans = exp(-alpha*(u-t))*r_t + m*(1-exp(-alpha*(u-t)))
    return ans
end

#int_S_T E_Q[r(u)|F_t]du:
integral_E_Q_r, _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), S,T)

integral_E_Q_r
#--------------------------------------------------------------------------#
# Calculating a_hat_3M:
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


function f_3M(t,S,T)
    " 
    Vasicek representation of f^{3M}(t,S,T) as described in Eq. (5.4) p.60
    "
    ans = (1/(T-S))*(exp(A(t,S,T) + B(t,S,T)*r_t) - 1)
    return ans
end

f_3M(0,S,T)

a_hat = integral_E_Q_r/((T-S)*f_3M(0,S,T))

#-------------------------------------------------------------------------------------------#
# 3M-arithmetic vs (a,b,c) 1M-SOFR futures: 
integral1 , _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), S,T1M)
integral2 , _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), T1M,T2M)
integral3 , _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), T2M,T)

f_1M_S_T1M = (1/(T1M-S))*integral1
f_1M_T1M_T2M = (1/(T2M-T1M))*integral2
f_1M_T2M_T = (1/(T-T2M))*integral3

#variable naming to be more consistent with MSc Thesis:
a = f_1M_S_T1M        #alpha, I use alpha in Vasicek, hence a: 
beta = f_1M_T1M_T2M   #beta
gamma = f_1M_T2M_T    #gamma

futures = [a,beta,gamma]

#E_Q[X^(3M_A)(S,T)|F_t] = q: 
q = (1/(T-S))*integral_E_Q_r 

#matrix of coeff: 
M = [a^(2)      a*beta  a*gamma;
     beta^(2)   a*beta  beta*gamma; 
     gamma^(2)  a*beta  beta*gamma] 
     
#vector of values: 
b = q.*[a;
        beta;
        gamma] 

#optimal weight of futures:
x_hat = inv(M)*b 
#-------------------------------------------------------------------------------------------#
# incase M is not invertible: 
# Define optimization problem
model = Model(HiGHS.Optimizer)
@variable(model, x[1:3])
@objective(model, Min, sum(x))
@constraint(model, M * x .== b)

# Solve optimization problem
optimize!(model)
# optimal value
x_tilde = value.(x)

#-------------------------------------------------------------------------------------------#
# Simulations: 
n_sim = 10^(6)
#constants: 
#int_S_T E_Q[r(u)|F_t]du:
integral_E_Q_r, _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), S,T)

futures_weighted_M_inv = x_hat'futures
futures_weighted_BP = x_tilde'futures


Random.seed!(1234)
X_3MA = zeros(n_sim)
for i in 1:n_sim
    #aritmetic interest rate relaization:
    X_3MA[i] = (1/(T-S))*(int_r_start_stop(S,T,t))
end

mean(X_3MA)
#elementwise substraction:
ER_1 = X_3MA .-(1/(T-S))*integral_E_Q_r
ER_2_M_inv = X_3MA .-futures_weighted_M_inv
ER_2_random = X_3MA .-[0.33, -0.33, 0.33]'futures

#-----------------------------------------------------------
# plotting of histograms: 
#hedge with a_{t}^{3M}- f^{3M}
mean_ER_1 = round(mean(ER_1), digits = 3)
sigma_ER_1 = round(std(ER_1), digits = 2) 

#hedge with optimal (a_{t}^{1M}, b_{t}^{1M}, c_{t}^{1M})-f^{1M}
mean_ER_2_M_inv = round(mean(ER_2_M_inv), digits = 3)
sigma_ER_2 = round(std(ER_2_M_inv), digits = 2)

#not optimal 1M hedges, naive strategy:  
mean_ER_2_random = round(mean(ER_2_random), digits = 3)
sigma_ER_2_random = round(std(ER_2_random), digits = 2)

#ER_1
histogram(ER_1, 
          color =:lightblue, 
          xlabel="Value", 
          ylabel="Frequency", 
          title = L"Histogram\; of\; ER_{1}(0), \; s_{ER_{1}}\approx 0.01,\;n_{sim} = 10^{6}", 
          labels = "ER_1(0)", 
          xticks = [-3*sigma_ER_1, -2*sigma_ER_1,-sigma_ER_1, mean_ER_1,  sigma_ER_1, 2*sigma_ER_1, 3*sigma_ER_1]
          )
vline!([mean_ER_1], lw = 5, labels = L"mean(ER_{1}(0))" )


#ER_2_M_inv:
histogram(ER_2_M_inv, 
          color =:lightblue, 
          xlabel="Value", 
          ylabel="Frequency", 
          title = L"Histogram\; of\; ER_{2}^{M_{inv}}(0), \; s_{ER_{2}^{M_{inv}}} \approx 0.01, \;n_{sim} = 10^{6}", 
          labels = L"ER_{2}^{M_{inv}}(0)", 
          xticks = [-3*sigma_ER_2, -2*sigma_ER_2,-sigma_ER_2, mean_ER_2_M_inv, sigma_ER_2, 2*sigma_ER_2, 3*sigma_ER_2]
          )
vline!([mean_ER_2_M_inv], lw = 5, labels = L"mean(ER_{2}^{M_{inv}}(0))")


#Naive strategy ER_2_random (0.33, -0.33, 0.33)
ticks = round.([mean_ER_2_random + i*sigma_ER_2_random for i in -2:1:2], digits = 3)

histogram(ER_2_random, 
          color =:lightblue, 
          xlabel="Value", 
          ylabel="Frequency", 
          title = L"Hist\; of\; ER_{2}^{(\hat{a}_{0}, \hat{b}_{0}, \hat{c}_{0})}(0), \; s_{ER_{2}^{(\hat{a}_{0}, \hat{b}_{0}, \hat{c}_{0})}} \approx 0.01, \;n_{sim} = 10^{6}", 
          labels = L"ER_{2}^{(0.33,-0.33,0.33)}(0)", 
          xticks = ticks
          )
vline!([mean_ER_2_random], lw = 5, labels =  L"mean(ER_{2}^{(0,0,1)}(0))")



#adressing normality:
#-----------------------------------------
x = ER_1[1:10^(6)]
y = rand(Normal(mean_ER_1, sigma_ER_1), 10^(6))

qqplot(x,y, title =  L"(Q-Q)\; plot\; of\; ER_{1}(0) \;vs\; \mathcal{N}\left(\overline{ER_1(0)}, s_{ER_1(0)}^{2})\right)", 
            xlabel = "Theoretical Quantiles", 
            ylabel = "Sample Quantiles")














