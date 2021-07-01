#==============================================================================#
# Asymmetric polymer/salt model system
#==============================================================================#

"""
    TwoChainModel{TC <: AbstractChainStructure} = 
        Union{AsymmetricCoacervate{TC}, AssociationCoacervate{TC}, SelfComplimentaryCoacervate{TC}}

Union type for models representing a standard coacervate model with treatment of A, C, (+), and (-) components.

Binodal variable format: [A1, C1, (+)1, (-)1, A2, C2, (+)2, (-)2, ν, Ψ]
"""
const TwoChainModel{TC <: AbstractChainStructure} = 
    Union{AsymmetricCoacervate{TC}, AssociationCoacervate{TC}, SelfComplimentaryCoacervate{TC}}

function neutralbulk(phi, model::TwoChainModel)
    sigA0, sigC0 = model.sig
    wA, wC, wP, wM = model.omega
    wS = (wP+wM)/2.0

    phiPS = phi[3] * wP / wS
    phiMS = phi[3] * wM / wS
    phiPC = sigA0 * phi[1] * wP / (wA*z)
    phiMC = sigC0 * phi[2] * wM / (wC*z)
    
    return [phi[1], phi[2], phiPS + phiPC, phiMS + phiMC]
end

function differencebulk(phi, model)
    phiA, phiC, phiS = phi
    wA, wC, wP, wM = model.omega
    sigA, sigC = model.sig

    chargeA = sigA*phiA/wA
    chargeC = sigC*phiC/wC
    diff = abs(chargeA - chargeC)/z

    if chargeA > chargeC
        return [phiA, phiC, phiS + wP*diff, phiS]
    else
        return [phiA, phiC, phiS, phiS + wM*diff]
    end
end

function chainstructs(model::TwoChainModel{TC}) where TC
    achain = ChainStructure{TC, eltype(model.dp)}(model.omega[1], model.dp[1], model.b[1], model.lp[1])
    cchain = ChainStructure{TC, eltype(model.dp)}(model.omega[2], model.dp[2], model.b[2], model.lp[2])
    return (achain, cchain)
end

function chainstructs(model::TwoChainModel{AdaptiveChain}, vars)
    achain = ChainStructure{AdaptiveChain, eltype(vars)}(model.omega[1], model.dp[1], model.b[1], vars[end-1])
    cchain = ChainStructure{AdaptiveChain, eltype(vars)}(model.omega[2], model.dp[2], model.b[2], vars[end])
    return (achain, cchain)
end

#==============================================================================#

bndlminx(state::BinodalState, model::TwoChainModel) = [state.dense[1], state.dense[2], 1-sum(state.dense), state.nu]

bndlsolvex(state::BinodalState, model::TwoChainModel) = [state.sup..., state.dense..., state.nu, get(state.props, :psi, 0.0)]

function bndlscale!(x, model::TwoChainModel)
    if length(x) == 4
        logscale!(x)
    elseif length(x) == 10
        logscale!(x, 1:9)
    end
    return nothing
end

function bndlunscale!(x, model::TwoChainModel)
    if length(x) == 4
        logunscale!(x)
    elseif length(x) == 10
        logunscale!(x, 1:9)
    end
    return nothing
end


function bndlstate(x, model::TwoChainModel)
    @assert length(x) == 4 || length(x) == 10 "Invalid number of parameters for binodal state. Requires (4, 10) for $(typeof(model))."

    if length(x) == 4
        phiAC, phiCC, phiWC, nu = x
        phiAB, phiCB, phiPB, phiMB = model.bulk
        phiWB = 1 - sum(model.bulk)
        wA, wC, wP, wM = model.omega
        sigA, sigC = model.sig
        
        # Derive dense phase params
        phiPC = wP/(wP+wM) - wP*(phiAC + phiCC + phiWC)/(wP+wM) + sigA*phiAC*(wM*wP/(z*wA*(wP+wM))) - sigC*phiCC*(wM*wP/(z*wC*(wP+wM)))
        phiMC = 1 - phiAC - phiCC - phiPC - phiWC
        
        # Derive supernatant phase params
        phiAS = (phiAB - nu*phiAC)/(1 - nu)
        phiCS = (phiCB - nu*phiCC)/(1 - nu)
        phiWS = (phiWB - nu*phiWC)/(1 - nu)
        phiPS = wP/(wP+wM) - wP*(phiAS + phiCS + phiWS)/(wP+wM) + sigA*phiAS*(wM*wP/(z*wA*(wP+wM))) - sigC*phiCS*(wM*wP/(z*wC*(wP+wM)))
        phiMS = 1 - phiAS - phiCS - phiPS - phiWS

        # Update state struct
        bulk = @SVector [phiAB, phiCB, phiPB, phiMB]
        sup = @SVector [phiAS, phiCS, phiPS, phiMS]
        dense = @SVector [phiAC, phiCC, phiPC, phiMC]
        
        return BinodalState(bulk, sup, dense, nu)
    else
        phiAB, phiCB, phiPB, phiMB = model.bulk
        phiAS, phiCS, phiPS, phiMS = x[1], x[2], x[3], x[4]
        phiAC, phiCC, phiPC, phiMC = x[5], x[6], x[7], x[8]
        nu, psi = x[9], x[10]

        bulk = @SVector [phiAB, phiCB, phiPB, phiMB]
        sup = @SVector [phiAS, phiCS, phiPS, phiMS]
        dense = @SVector [phiAC, phiCC, phiPC, phiMC]

        state = BinodalState(bulk, sup, dense, nu)
        state.props[:psiG] = x[10]

        return state
    end
end

function bndlf!(F, xs, model::TwoChainModel; scaled::Bool = false)
    # Scale coordinates
    x = scaled ? bndlunscale(xs, model) : xs

    # Offload parameters
    phiAB, phiCB, phiPB, phiMB = model.bulk
    phiAS, phiCS, phiPS, phiMS = x[1], x[2], x[3], x[4]
    phiAC, phiCC, phiPC, phiMC = x[5], x[6], x[7], x[8]
    nu, psi = x[9], x[10]
    
    sup = @SVector [phiAS, phiCS, phiPS, phiMS]
    coac = @SVector [phiAC, phiCC, phiPC, phiMC]
    wA, wC, wP, wM = model.omega
    sigA, sigC = model.sig
    
    muS, pS = mupressure(sup, model)
    muC, pC = mupressure(coac, model)

    F[1] = muS[1] - muC[1] + sigA*psi/wA
    F[2] = muS[2] - muC[2] - sigC*psi/wC
    F[3] = muS[3] - muC[3] - z*psi/wP
    F[4] = muS[4] - muC[4] + z*psi/wM
    F[5] = pS - pC
    F[6] = sigC*phiCS/wC - sigA*phiAS/wA + z*phiPS/wP - z*phiMS/wM
    F[7] = sigC*phiCC/wC - sigA*phiAC/wA + z*phiPC/wP - z*phiMC/wM
    F[8] = phiAB - (1 - nu)*phiAS - nu*phiAC
    F[9] = phiCB - (1 - nu)*phiCS - nu*phiCC
    F[10] = phiPB - (1 - nu)*phiPS - nu*phiPC

    return nothing
end

function bndlj!(J, xs, model::TwoChainModel; scaled::Bool = false)
    # Scale coordinates
    x = scaled ? bndlunscale(xs, model) : xs

    # Offload parameters
    phiAB, phiCB, phiPB, phiMB = model.bulk
    phiAS, phiCS, phiPS, phiMS = x[1], x[2], x[3], x[4]
    phiAC, phiCC, phiPC, phiMC = x[5], x[6], x[7], x[8]
    nu, psi = x[9], x[10]
    
    sup = [phiAS, phiCS, phiPS, phiMS]
    coac = [phiAC, phiCC, phiPC, phiMC]
    sig = [-model.sig[1], model.sig[2], z, -z]
    omega = model.omega

    # Fill the Jacobian, J = ∂F(x)/∂x
    # Much of Jacobian information is contained Hessian of free energy (∂mu/∂phi terms == 2nd deriv. free energy)
    # Remaining terms are simple to extract from F(x) definition
    fill!(J, zero(eltype(J)))

    f2S = f2total(sup, model)
    f2C = f2total(coac, model)

    # 1st 4 rows of J == ∂mu/∂phi (negative for coacervate)
    # Also constant term from psi derivative in 10th column
    for i = 1:4
        for k = 1:4
            J[i,k] = f2S[i,k]
            J[i,k+4] = -f2C[i,k]
        end
        J[i,10] = -sig[i]/omega[i]
    end

    # 5th row J == ∂P/∂phi --> J[3,k] = dot(phi_i ∂mu_i/∂phi_k)
    for k = 1:4
        J[5,k] = dot(sup, f2S[:,k])
        J[5,k+4] = -dot(coac, f2C[:,k])
    end

    # 6th/7th rows from electroneutrality
    for k = 1:4
        J[6,k] = sig[k]/omega[k]
        J[7,k+4] = sig[k]/omega[k]
    end

    # Jacobian 8/9/10 rows from mass balance terms
    # 8  = {nu - 1, 0, 0, 0, -nu, 0, 0, 0, phiAS - phiAC, 0}
    J[8,1] = nu - 1
    J[8,5] = -nu
    J[8,9] = phiAS - phiAC

    # 9  = {0, nu - 1, 0, 0, 0, -nu, 0, 0, phiCS - phiCC, 0}
    J[9,2] = nu - 1
    J[9,6] = -nu
    J[9,9] = phiCS - phiCC

    # 10 = {0, 0, nu - 1, 0, 0, 0, -nu, 0, phiPS - phiPC, 0}
    J[10,3] = nu - 1
    J[10,7] = -nu
    J[10,9] = phiPS - phiPC

    # Scale Jacobian by variable transformation chain rule, if needed
    if scaled
        dx = [exp(xs[i])/(1 + exp(xs[i]))^2 for i = 1:9]
        push!(dx, 1.0) # Galvani potential is never scaled, remains unity

        for i in 1:size(J,1)
            for k in 1:size(J,2)
                J[i,k] *= dx[k]
            end
        end
    end
    
    return nothing
end
