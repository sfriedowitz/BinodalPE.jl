"""
    AsymmetricCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}
    
A coacervate model with asymmetric treatment of oppositely charged polyions and small salt ions.
"""
mutable struct AsymmetricCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}
    bulk   :: MVector{4,Float64}
    omega  :: SVector{4,Float64}
    smear  :: SVector{4,Float64}
    sig    :: SVector{2,Float64}
    chi    :: SVector{3,Float64}
    np     :: SVector{2,Float64}
    b      :: SVector{2,Float64}
    lp     :: SVector{2,Float64} # Only used if chain is WormLike
    lB     :: Float64
    vargs  :: Dict{Symbol,Any}
end

function AsymmetricCoacervate(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0, 1.0, 1.0])
    sig = get(kwargs, :sig, [1.0, 1.0])
    chi = get(kwargs, :chi, [0.0, 0.0, 0.0])
    np = get(kwargs, :np, [100.0, 100.0])
    b = get(kwargs, :b, [1.0, 1.0])
    lp = get(kwargs, :lp, [1.0, 1.0])
    lB = get(kwargs, :lB, lBbar)
    vargs = get(kwargs, :vargs, Dict())

    smear = 0.5 .* omega .^ (1/3)
    model = AsymmetricCoacervate{structure}(zeros(4), omega, smear, sig, chi, np, b, lp, lB, vargs)
    return model
end

#==============================================================================#

function free_energy(phi, model::AsymmetricCoacervate)
    return ftranslational(phi, model) + fchi(phi, model) + felectrostatic(phi, model)
end

function free_energy(phi, model::AsymmetricCoacervate{AdaptiveChain})
    vars = varsolve(phi, model; model.vargs...)
    return ftranslational(phi, model) + fchi(phi, model) + felectrostatic(phi, vars, model)
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

varscale(x, model::AsymmetricCoacervate{AdaptiveChain}) = log.(x)

varunscale(x, model::AsymmetricCoacervate{AdaptiveChain}) = exp.(x)

varinit(phi::AbstractArray{TF}, model::AsymmetricCoacervate{AdaptiveChain}) where TF = TF.(Vector(model.b/2))

function varf!(F::AbstractVector{TF}, x, phi, model::AsymmetricCoacervate{AdaptiveChain}) where TF
    vars = varunscale(x, model)
    F[1], F[2] = dftot_adaptive(phi, vars, model)
    return nothing
end