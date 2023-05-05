using Revise

#stats etc.
using Random 
using Distributions

#plotting histogram and LaTeX labels:
using Plots 
using LaTeXStrings

#Vasicek dynamics of SOFR futures rates:
alpha = 0.30
sigma = 0.03

function B(t,S,T)
    return (1/alpha)*(exp(-alpha*(S-t))-exp(-alpha*(T-t)))
end


function f1M(t,S,T, dt, initial)
    "
    Simulates 1M-futures rates for t<= S
    Args:
        t {Float64}: initial start point 
        S {Float64}: start observation period 
        T {Float65}: end observation period 
        dt {Float64}: stepsize 
        inital {Float64}: inital futures rate
    Returns: 
        df1M(t,S,T) = 1/(T-S)*B(t,S,T)*sigma*dW^{Q}(t)
    "
    time = range(t, S, step=dt)
    n = length(time)

    Random.seed!(1)
    W = rand(Normal(0,1), n)
    f1_r = zeros(n)
    f1_r[1] = initial
    for i in 2:n
        f1_r[i] = f1_r[i-1] + (1/(T-S))*B(time[i-1], S,T)*sigma*dt*W[i]
    end

    return f1_r
end


function f3M(t,S,T, dt, initial)
    "
    Simulates 3M-futures rates for t<= S
    Args:
        t {Float64}: initial start point 
        S {Float64}: start observation period 
        T {Float65}: end observation period 
        dt {Float64}: stepsize 
        inital {Float64}: inital futures rate
    Returns: 
        df3M(t,S,T) = (f3M(t,S,T)+1/(T-S))*B(t,S,T)*sigma*dW^{Q}(t)
    "
    time = range(t, S, step=dt)
    n = length(time)

    Random.seed!(2)
    W = rand(Normal(0,1), n)
    f3_r = zeros(n)
    f3_r[1] = initial
    for i in 2:n
        f3_r[i] = f3_r[i-1] + (f3_r[i-1] + 1/(T-S))*B(time[i-1], S,T)*sigma*dt*W[i]
    end

    return f3_r
end

Random.seed!(3)
#time params: 
t = 0
dt = 1/360
#1M-futures:
S1M = 6/12
T1M = S1M + 1/12
#3M-futures: 
S3M = S1M
T3M = S3M + 3/12


#simulation 
time = range(t, S1M, step=dt)
n = length(time)
W = rand(Normal(0,1), n)
f1 = zeros(n)
f3 = zeros(n)


f1[1] = (100-95.025)*1/100
f3[1] = (100-95.16)*1/100

for i in 2:n
    f1[i] = f1[i-1] + (1/(T1M-S1M))*B(time[i-1], S1M,T1M)*sigma*dt*W[i]
    f3[i] = f3[i-1] + (f3[i-1] + 1/(T3M-S3M))*B(time[i-1], S3M,T3M)*sigma*dt*W[i]
end

plot(f1, label = L"f^{1M}(t, S_{1M},T_{1M})", title = L"\alpha = 0.30,\; \sigma = 0.03,\; t\in [0,S_{1M}]", legend= :topleft)
plot!(f3, label = L"f^{3M}(t, S_{3M}, T_{3M})")
xticks!([0, n/2 ,n], ["0", L"\frac{S_{1M}}{2}", L"S_{1M}"])


#=
f1M_initial = (100-95.025)*1/100
f3M_initial = (100-95.16)*1/100

time = range(t, S1M, step=dt)
n = length(time)


plot(f1M(0, S1M, T1M, dt,f1M_initial), label = L"f^{1M}(t, S_{1M},T_{1M})", title = L"\alpha = 0.30,\; \sigma = 0.03,\; t\in [0,S_{1M}]", legend= :topleft)
plot!(f3M(0, S3M, T3M, dt, f3M_initial), label = L"f^{3M}(t, S_{3M}, T_{3M})")
xticks!([0, n/2 ,n], ["0", L"\frac{S_{1M}}{2}", L"S_{1M}"])
=#
