using GREMLModels
using SnpArrays, DataFrames, CSV, LinearAlgebra
using StatsBase, StatsModels
using JLD
using Distributions

idges = Symmetric(idges, :U)
iges = Symmetric(iges, :U)
dges = Symmetric(dges, :U)
S = Symmetric(S, :U)

y = yX[:,1]
X = ones(length(y))

# FULL MODEL

r = [dges, iges, idges, S, Diagonal(ones(length(y)))]
 
dat = GREMLData(y, X, r)

lbs = [0.0, -Inf, 0.0, 0.0, 0.0]
ini =  sqrt.(var(y) .* [0.3, 0.0, 0.0, 0.0, 0.7])
@time m = GREMLModel(dat, ini, lbs, false)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2 + θ[3]^2
    δ[3] = θ[1] * θ[2]
    δ[4] = θ[4]
    δ[5] = θ[5]
    δ
end
@time fit!(m)

m.δ

hessian!(m.opt.H, m)
j = jacobian(m)
vcovvc(m)
vc_delta = j * vcovvc(m) * j'
sqrt.(diag(vc_delta))

# SEMI-REDUCED MODEL
rsemi = [dges, iges, S, Diagonal(ones(length(y)))]
 
function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]^2
    δ[3] = θ[3]
    δ[4] = θ[4]
    δ
end

datsemi = GREMLData(y, X, rsemi)

lbs = [0.0, 0.0, 0.0, 0.0]
ini =  sqrt.(var(y) .* [0.3, 0.0, 0.0, 0.7])
@time msemi = GREMLModel(datsemi, ini, lbs, false)

@time fit!(msemi)

msemi.δ

hessian!(msemi.opt.H, msemi)
jsemi = jacobian(msemi)
vcovvc(msemi)
vc_delta_semi = jsemi * vcovvc(msemi) * jsemi'
sqrt.(diag(vc_delta_semi))

# REDUCED MODEL

rred = [dges, S, Diagonal(ones(length(y)))]
 
datred = GREMLData(y, X, rred)
function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    δ[3] = θ[3]
    δ
end
lbs = [0.0, 0.0, 0.0]
ini =  sqrt.(var(y) .* [0.3, 0.0, 0.7])
@time mred = GREMLModel(datred, ini, lbs, false)

function GREMLModels.transform!(δ::Vector, θ::Vector)
    δ[1] = θ[1]^2
    δ[2] = θ[2]
    δ[3] = θ[3]
    δ
end
@time fit!(mred)

mred.δ

hessian!(mred.opt.H, mred)
jred = jacobian(mred)
vcovvc(mred)
vc_delta_red = jred * vcovvc(mred) * jred'
sqrt.(diag(vc_delta_red))

# RESULTS

m.δ
sqrt.(diag(vc_delta))

msemi.δ
sqrt.(diag(vc_delta_semi))

mred.δ
sqrt.(diag(vc_delta_red))

chi_m_msemi = msem.opt.ffinal-m.opt.ffinal
chi_msemi_mred =  mred.opt.ffinal-msemi.opt.ffinal

chi_dist = Chisq(1)
p_m_msemi = ccdf(chi_dist, chi_m_msemi)
p_msemi_mred = ccdf(chi_dist, chi_msemi_mred)
