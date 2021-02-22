"""
    SelfComplimentaryCoacervate{TC <: AbstractChainStucture} <: AbstractModel{TC}

An asymmetric coacervate system with treatment of reversible ionic association
and non-ionic, self-complimentary association between polyanions.
The self-complimentary binding is intended to model pairing of DNA base regions.
"""
mutable struct SelfComplimentaryCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}
    bulk   :: MVector{4,Float64}
    omega  :: SVector{4,Float64}
    sig    :: SVector{2,Float64}
    fA     :: Float64 # Base-fraction of self-complimentary stickers on polyanion
    vAA    :: Float64 # Reduced volume of an entire self-complimentary bond 
    dgAA   :: Float64
    assoc  :: AssociationCoacervate{TC}
end

function SelfComplimentaryCoacervate(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0, 1.0, 1.0])
    sig = get(kwargs, :sig, [1.0, 1.0])
    fA = get(kwargs, :fA, 1.0)
    vAA = get(kwargs, :vAA, 1.0)
    dgAA = get(kwargs, :dgAA, 0.0)
    dg = get(kwargs, :dg, [0.0, 0.0, 0.0])

    # Use association model struct as field
    assoc = AssociationCoacervate(; structure = structure, dg = dg, kwargs...)
    model = SelfComplimentaryCoacervate{structure}(zeros(4), omega, sig, fA, vAA, dgAA, assoc)

    return model
end

#==============================================================================#

function ftotal(phi, model::SelfComplimentaryCoacervate)
    # Solve for all assoc vars (including self-comp term)
    theta = selfcompsolve(phi, model)
    assoc = varinit(phi, model.assoc)
    try
        assoc = varsolve(phi, model.assoc)
    catch e
        @error "Association solver failure." maxlog = 1
        nothing
    end

    # Majority of free energy from association sub-model
    ftrans = ftranslational(phi, assoc, model.assoc)
    fcomb = fcombinatorial(phi, assoc, model.assoc)
    finf = fbinding(phi, assoc, model.assoc)
    ffh = fchi(phi, model.assoc)
    fel = felectrostatic(phi, assoc, model.assoc)

    # Self-comp part of free energy
    fself = fselfcomp(phi, theta, model)

    return ftrans + fcomb + finf + ffh + fel + fself
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

function selfcompsolve(phi, model::SelfComplimentaryCoacervate)
    # Analytically solution to the self-complimentary fraction
    phiA = phi[1]
    wA = model.omega[1]
    fA = model.fA
    vAA = model.vAA
    dgAA = model.dgAA
    kAA = exp(-dgAA)

    cA = fA*phiA*vAA*kAA/wA
    if cA < 1e-8
        tA = cA
    else
        tA = 1.0 + 1/(2.0*cA) - sqrt(1.0 + 4.0*cA)/(2.0*cA)
    end

    return tA
end

function varsolve(phi, model::SelfComplimentaryCoacervate)
    # Solve for both binding parts
    vars = varsolve(phi, model.assoc)
    tA = selfcompsolve(phi, model)
    push!(vars, tA)
    return vars
end