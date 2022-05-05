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
    chi    :: Float64
    dp     :: Float64
    b      :: Float64
    lp     :: Float64
    lB     :: Float64
    zP     :: Float64
    zM     :: Float64
    vargs  :: Dict{Symbol,Any}
end

function SinglePolyion(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0, 1.0])
    sig = get(kwargs, :sig, 1.0)
    dg = get(kwargs, :dg, 0.0)
    chi = get(kwargs, :chi, 0.0)
    dp = get(kwargs, :dp, 100)
    b = get(kwargs, :b, 1.0)
    lp = get(kwargs, :lp, 1.0)
    lB = get(kwargs, :lB, lBbar)
    zP = get(kwargs, :zP, 1.0)
    zM = get(kwargs, :zM, 1.0)
    vargs = get(kwargs, :vargs, Dict())

    smear = 0.5 .* omega .^ (1/3)
    model = SinglePolyion{structure}(zeros(3), omega, smear, sig, dg, chi, dp, b, lp, lB, zP, zM, vargs)
    return model
end

#==============================================================================#

chainstructs(model::SinglePolyion{TC}) where TC = ChainStructure{TC, typeof(model.dp)}(model.omega[1], model.dp, model.b, model.lp)

chainstructs(model::SinglePolyion{AdaptiveChain}, vars) = ChainStructure{AdaptiveChain, eltype(vars)}(model.omega[1], model.dp, model.b, vars[2])

function neutralbulk(phi, model::SinglePolyion)
    @unpack sig, omega, zP = model
    wA, wP, wM = omega
    
    counter = sig * phi[1] * wP / (wA*zP)
    return [phi[1], phi[2] + counter, phi[2]]
end

function ftotal(phi, model::SinglePolyion)
    assoc = varsolve(phi, model; model.vargs...)

    ftot = ftranslational(phi, assoc, model)
    ftot += fchi(phi, model)
    ftot += fcombinatorial(phi, assoc, model)
    ftot += fbinding(phi, assoc, model)
    ftot += felectrostatic(phi, assoc, model)

    return ftot
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

function varscale(x, model::SinglePolyion{AdaptiveChain})
    new = similar(x)
    new[1] = logscale(x[1])
    new[2] = log(x[2])
    return new
end

function varunscale(x, model::SinglePolyion{AdaptiveChain})
    new = similar(x)
    new[1] = logunscale(x[1])
    new[2] = exp(x[2])
    return new
end

function varinit(phi, model::SinglePolyion)
    phiA, phiP, phiM = phi
    @unpack dg, b, zP = model

    # effective equilibrium constant
    muel = muel_association(phi, [0.5, 2.5], model)
    kA = exp(-dg - zP*muel + 1)

    # lumped constant for solution
    inner = (kA*phiP)^(1/zP)
    
    alpha0 = inner < Inf ? inner/(1 + inner) : 1.0
    
    alpha0 = clamp(alpha0, 1e-10, 0.99)
    
    init = [alpha0]
    if isa(model, SinglePolyion{AdaptiveChain})
        push!(init, b/2) # This just works better tbh
    end

    return init
end


#==============================================================================#

function varf!(F::AbstractVector{TF}, x, phi, model::SinglePolyion) where TF
    phiA, phiP, phiM = phi
    @unpack omega, dg, zP = model
    wA, wP, wM = omega

    alpha = exp(x[1])/ (1 + exp(x[1]))
    phiPF = phiP - alpha*phiA*wP/(wA*zP)
    sig = 1 - alpha

    muel = muel_association(phi, [alpha], model)
    F[1] = log((alpha/sig)^zP/phiPF) + dg + zP*muel - 1

    return nothing
end

function varf!(F::AbstractVector{TF}, x, phi, model::SinglePolyion{AdaptiveChain}) where TF
    phiA, phiP, phiM = phi
    @unpack omega, dg, zP = model
    wA, wP, wM = omega

    vars = varunscale(x, model)
    alpha, _ = vars

    phiPF = phiP - alpha*phiA*wP/(wA*zP)
    sig = 1 - alpha

    muel = muel_association(phi, vars, model)
    dftot = dftot_adaptive(phi, vars, model)

    F[1] = log((alpha/sig)^zP/phiPF) + dg + zP*muel - 1
    F[2] = dftot

    return nothing
end