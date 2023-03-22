using Random 
using Distributions
using Statistics

using LinearAlgebra
using Queryverse
using DataFrames

#for numerical integration:
using QuadGK

#plotting histogram:
using Plots 

#Vasicek parameters: 
alpha = 0.25
m = 0.035
r_t = 0.05
sigma = 0.1

#time parameters
t = 0
S = 1/12
T = 4/12
n_steps = 10000

function r_Vasicek(alpha, m,sigma, r_t, time_interval, n_steps)
    """
    Args: 
        #Vasicek parameters:
        alpha (Float): speed of reversion 
        m (float): long term mean level 
        sigma (float): volatility 
        r_t (float): initial value of r = (r(u))
        
        #time:
        time_interval (vector): the time interval we model over e.g [0,3]
        n_steps (int): number of timesteps we partition over
    
    Returns: 
        it simulates the process: r = (r(u)), for t in [t_start, t_end]
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
        r[i] = r[i-1] - alpha*(m-r[i-1])*dt + sigma*Z[i]*dt
    end

    return r
end


function int_r_start_stop(lower, upper, t)
    """
    Args: 
        lower (Float): where the integral starts 
        uppper (Float): where the integral stops
        t (float): when the simulation of r should start,
                   effects r_lower. 
    Returns: 
        the integral int_lower_upper r(u)du
    """

    function B(u)
        #=
        actually a function B(upper-u), of two variables, 
        we parametrize by one variable for easing integral in Julia
        =#
            return (1/alpha)*(1-exp(-alpha*(upper-u)))
    end

    function b(u)
        # again actually b(upper-u). 
        ans = -(m/alpha)*(1- exp(-alpha*(upper-u))) - m*(upper-u)
        return ans
    end

    #represents the integrand in N(0, int_lower_upper c^2(u)du)
    function c_2(u)
        integrand = (sigma/alpha)^(2)*(exp(-alpha*(upper-u)-1))^(2)
        return integrand
    end

    r_lower = last(r_Vasicek(alpha, m, sigma,r_t, [t,lower], n_steps))

    integral_c_2, _ = quadgk(c_2, lower,upper)

    ans = B(lower)*r_lower + b(lower) - (integral_c_2^(1/2))*rand(Normal(0,1))

    return ans
end 

int_r_start_stop(S, T, t)

function integrand_E_Q_r(u, r_t, t)
    ans = exp(-alpha*(u-t))*r_t + m*(1-exp(-alpha*(u-t)))
    return ans
end

#int_S_T E_Q[r(u)|F_t]du:
integral_E_Q_r, _ = quadgk(u -> integrand_E_Q_r(u, r_t, t), S,T)
quadgk(u -> integrand_E_Q_r(u, r_t, t), S,T)[1]

integral_E_Q_r

n_sim = 10^(6)
Er = zeros(n_sim)
for i in 1:n_sim
    Er[i] = (1/(T-S))*(int_r_start_stop(S,T,t) - integral_E_Q_r)
end


mean_Er = mean(Er)
mean_Er = round(mean_Er; digits = 3)
std_Er = std(Er)

histogram(Er, 
          color =:lightblue, 
          xlabel="Value", 
          ylabel="Frequency", 
          title="Histogram of Er(0), n_sim = 10^(6)", 
          labels = "Er(0)")
xlabel!("Er(0)")
plot!([mean_Er], seriestype="vline",
      xticks =([-1.5, -1.0, 0.5, mean_Er,0.5, 1.0, 1.5]), 
      label="mean(Er(0))", lw = 5, 
      dpi = 350)


#save plot
savefig("C:\\Users\\Andre\\OneDrive\\Dokumenter\\UiO\\Master fag\\Masteroppgave Fred\\Masteroppgave_Julia\\Error_hist.png")

round(mean_Er; digits = 3)


histogram(abs.(Er))
