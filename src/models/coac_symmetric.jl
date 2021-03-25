"""
    SymmetricCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}

A coacervate model with a symmetric treatment of oppositely charged polyions and small salt ions.
"""
mutable struct SymmetricCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}
    bulk   :: MVector{2,Float64}
    omega  :: SVector{2,Float64}
    smear  :: SVector{2,Float64}
    sig    :: Float64
    chi    :: Float64
    dp     :: Float64
    b      :: Float64
    lp     :: Float64
    lB     :: Float64
    vargs  :: Dict{Symbol,Any}
end

function SymmetricCoacervate(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0])
    sig = get(kwargs, :sig, 1.0)
    chi = get(kwargs, :chi, 0.0)
    dp = get(kwargs, :dp, 100)
    b = get(kwargs, :b, 1.0)
    lp = get(kwargs, :lp, 1.0)
    lB = get(kwargs, :lB, lBbar)
    vargs = get(kwargs, :vargs, Dict())

    smear = 0.5 .* omega .^ (1/3)
    model = SymmetricCoacervate{structure}(zeros(2), omega, smear, sig, chi, dp, b, lp, lB, vargs)
    return model
end

#==============================================================================#

chainstructs(model::SymmetricCoacervate{TC}) where TC = ChainStructure{TC}(model.dp, model.lp, model.b, model.omega[1])

chainstructs(model::SymmetricCoacervate{AdaptiveChain}, vars) = ChainStructure{AdaptiveChain}(model.dp, vars[1], model.b, model.omega[1])

ftotal(phi, model::SymmetricCoacervate) = ftranslational(phi, model) + fchi(phi, model) + felectrostatic(phi, model)

function ftotal(phi, model::SymmetricCoacervate{AdaptiveChain})
    vars = varsolve(phi, model; model.vargs...)
    ftranslational(phi, model) + fchi(phi, model) + felectrostatic(phi, vars, model)
end

#==============================================================================#
# Binodal solver functions
#==============================================================================#

bndlminx(state::BinodalState, model::SymmetricCoacervate) = [state.dense..., state.nu]

bndlsolvex(state::BinodalState, model::SymmetricCoacervate) = [state.sup..., state.dense..., state.nu]

bndlscale!(x, model::SymmetricCoacervate) = logscale!(x)

bndlunscale!(x, model::SymmetricCoacervate) = logunscale!(x)

function bndlstate(x, model::SymmetricCoacervate)
    @assert length(x) == 3 || length(x) == 5 "Invalid number of parameters for binodal state. Requires (3, 5) for $(typeof(model))."

    if length(x) == 3
        phiPC, phiSC, nu = x
        phiPB, phiSB = model.bulk
        phiWB = 1 - sum(model.bulk)

        # Get supernatant parameters
        phiPS = (phiPB - nu*phiPC)/(1 - nu)
        phiSS = (phiSB - nu*phiSC)/(1 - nu)
    
        # Update state struct
        bulk = @SVector [phiPB, phiSB]
        sup = @SVector [phiPS, phiSS]
        dense = @SVector [phiPC, phiSC]

        return BinodalState(bulk, sup, dense, nu)
    else
        phiPB, phiSB = model.bulk
        phiPS, phiSS, phiPC, phiSC, nu = x

        bulk = @SVector [phiPB, phiSB]
        sup = @SVector [phiPS, phiSS]
        dense = @SVector [phiPC, phiSC]

        return BinodalState(bulk, sup, dense, nu)
    end
end

function bndlf!(F, xs, model::SymmetricCoacervate; scaled::Bool = false)
    # Unscale coordinates
    x = scaled ? bndlunscale(xs, model) : xs

    # Offload parameters
    phiPB, phiSB = model.bulk
    phiPS, phiSS, phiPC, phiSC, nu = x
    
    sup = @SVector [phiPS, phiSS]   
    dense = @SVector [phiPC, phiSC]
    muS, pS = mupressure(sup, model)
    muC, pC = mupressure(dense, model)

    F[1] = muS[1] - muC[1]
    F[2] = muS[2] - muC[2]
    F[3] = pS - pC
    F[4] = phiPB - nu*phiPC - (1-nu)*phiPS
    F[5] = phiSB - nu*phiSC - (1-nu)*phiSS
    
    return nothing
end

function bndlj!(J, xs, model::SymmetricCoacervate; scaled::Bool = false)
    # Unscale coordinates
    x = scaled ? bndlunscale(xs, model) : xs

    # Offload parameters
    phiPB, phiSB = model.bulk
    phiPS, phiSS, phiPC, phiSC, nu = x

    sup = @SVector [phiPS, phiSS]   
    coac = @SVector [phiPC, phiSC]

    # Fill the Jacobian
    fill!(J, zero(eltype(J)))
    
    f2S = f2total(sup, model)
    f2C = f2total(coac, model)

    # First two rows J == ∂mu/∂phi (negative for coacervate)
    for i = 1:2
        for k = 1:2
            J[i,k] = f2S[i,k]
            J[i,k+2] = -f2C[i,k]
        end
    end

    # 3rd row J == ∂P/∂phi --> J[3,k] = dot(phi_i ∂mu_i/∂phi_k)
    for k = 1:2
        J[3,k] = dot(sup, f2S[:,k])
        J[3,k+2] = -dot(coac, f2C[:,k])
    end

    # 4th row J == {nu-1, 0, -nu, 0, phiPS-phiPC}
    J[4,1] = nu - 1
    J[4,3] = -nu
    J[4,5] = phiPS - phiPC

    # 5th row J == {0, nu-1, 0, -nu, phiSS-phiSC}
    J[5,2] = nu - 1
    J[5,4] = -nu
    J[5,5] = phiSS - phiSC

    # Scale Jacobian by variable transformation chain rule, if needed
    if scaled
        dx = [exp(s)/(1 + exp(s))^2 for s in xs]
        for i in 1:size(J,1)
            for k in 1:size(J,2)
                J[i,k] *= dx[k]
            end
        end
    end

    return nothing
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

varscale(x, model::SymmetricCoacervate{AdaptiveChain}) = log.(x)

varunscale(x, model::SymmetricCoacervate{AdaptiveChain}) = exp.(x)

varinit(phi::AbstractArray{TF}, model::SymmetricCoacervate{AdaptiveChain}) where TF = TF.([model.b/2])

function varf!(F, x, phi, model::SymmetricCoacervate{AdaptiveChain})
    vars = varunscale(x, model)
    F[1] = dftot_adaptive(phi, vars, model)
    return nothing
end