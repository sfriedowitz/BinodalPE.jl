"""
    SinglePolyion{TC <: AbstractChainStructure} <: AbstractModel{TC}

A model for a single polyanion solution in the presence of like- and oppositely-charged small ions.

Oppositely charged small ions are capable of condensing on the polyion backbone.
The extent of such condensation is considered by the reversible binding reactions in `varsolve`.

The composition vector `phi` is ordered as: [A, (+), (-)]
"""
mutable struct SinglePolyion{TC <: AbstractChainStructure} <: AbstractModel{TC}
    bulk   :: MVector{3,Float64}
    omega  :: SVector{3,Float64}
    smear  :: SVector{3,Float64}
    sig    :: Float64
    dg     :: Float64
    np     :: Float64
    chi    :: Float64
    b      :: Float64
    lp     :: Float64
    lB     :: Float64
    vargs  :: Dict{Symbol,Any}
end

function SinglePolyion(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0, 1.0])
    sig = get(kwargs, :sig, 1.0)
    dg = get(kwargs, :dg, 0.0)
    np = get(kwargs, :np, 100)
    chi = get(kwargs, :chi, 0.0)
    b = get(kwargs, :b, 1.0)
    lp = get(kwargs, :lp, 1.0)
    lB = get(kwargs, :lB, lBbar)
    vargs = get(kwargs, :vargs, Dict())

    smear = 0.5 .* omega .^ (1/3)
    model = SinglePolyion{structure}(zeros(3), omega, smear, sig, dg, np, chi, b, lp, lB, vargs)
    return model
end

#==============================================================================#

function neutralbulk(phi, model::SinglePolyion)
    sig = model.sig
    wA, wP, wM = model.omega
    
    counter = sig * phi[1] * wP / wA
    return [phi[1], phi[2] + counter, phi[2]]
end

function ftranslational(phi, model::SinglePolyion)
    assoc = varsolve(phi, model)
    return ftranslational(phi, assoc, model)
end

function free_energy(phi, model::SinglePolyion)
    assoc = varsolve(phi, model; model.vargs...)

    ftot = ftranslational(phi, assoc, model)
    ftot += fchi(phi, model)
    ftot += fcombinatoric(phi, assoc, model)
    ftot += fbinding(phi, assoc, model)
    ftot += felectrostatic(phi, assoc, model)

    return ftot
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

function varscale(x, model::SinglePolyion{AdaptiveChain})
    new = similar(x)
    new[1] = log(x[1]/(1 - x[1]))
    new[2] = log(x[2])
    return new
end

function varunscale(x, model::SinglePolyion{AdaptiveChain})
    new = similar(x)
    new[1] = exp(x[1])/(1 + exp(x[1]))
    new[2] = exp(x[2])
    return new
end

function varinit(phi, model::SinglePolyion)
    phiA, phiP, phiM = phi
    dg = model.dg
    b = model.b

    muel = muel_assoc(phi, [0.5, 2.5], model)
    kA = exp(-dg - muel + 1)

    if phiP < 1e-10 # To avoid roundoff error
        alpha0 = phiP*kA - 2*(phiP*kA)^2
    else
        inner = 1 + kA^2*(phiA-phiP)^2 + 2*kA*(phiA+phiP)
        quad = inner < 0 ? 1 : sqrt(inner)
        alpha0 = 0.5 + 0.5/(kA*phiA) + phiP/(2*phiA) - quad/(2*kA*phiA)
    end
    alpha0 = clamp(alpha0, 1e-10, 0.95)
    
    init = [alpha0]
    if isa(model, SinglePolyion{AdaptiveChain})
        push!(init, b/2) # This just works better tbh
    end

    return init
end

#==============================================================================#

function varf!(F::AbstractVector{TF}, x, phi, model::SinglePolyion) where TF
    phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    dg = model.dg

    alpha = exp(x[1])/ (1 + exp(x[1]))
    phiPF = phiP - alpha*phiA*wP/wA
    sig = 1 - alpha

    muel = muel_assoc(phi, [alpha], model)
    F[1] = log(alpha/(1-alpha)/phiPF) + dg + muel - 1

    return nothing
end

function varf!(F::AbstractVector{TF}, x, phi, model::SinglePolyion{AdaptiveChain}) where TF
    phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    dg = model.dg

    vars = varunscale(x, model)
    alpha, lp = vars

    phiPF = phiP - alpha*phiA*wP/wA
    sig = 1 - alpha

    muel = muel_assoc(phi, vars, model)
    dftot = dftot_adaptive(phi, vars, model)

    F[1] = log(alpha/(1-alpha)/phiPF) + dg + muel - 1
    F[2] = dftot

    return nothing
end