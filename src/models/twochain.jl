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
    @unpack sig, omega, zP, zM = model
    sigA0, sigC0 = sig
    wA, wC, wP, wM = omega
    wS = (zM*wP+zP*wM)/(zP + zM)

    phiPS = phi[3] * wP / wS
    phiMS = phi[3] * wM / wS
    phiPC = sigA0 * phi[1] * wP / (wA*zP)
    phiMC = sigC0 * phi[2] * wM / (wC*zM)
    
    return [phi[1], phi[2], phiPS + phiPC, phiMS + phiMC]
end

function differencebulk(phi, model)
    phiA, phiC, phiS = phi
    @unpack omega, sig, zP, zM = model
    wA, wC, wP, wM = omega
    sigA, sigC = sig

    chargeA = sigA*phiA/wA
    chargeC = sigC*phiC/wC
    diff = abs(chargeA - chargeC)

    if chargeA > chargeC
        return [phiA, phiC, phiS + wP*diff/zP, phiS]
    else
        return [phiA, phiC, phiS, phiS + wM*diff/zM]
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

    x0 = similar(x)
    phiAS, phiCS, phiPS, phiMS, phiAC, phiCC, phiPC, phiMC, nu, psi = x
    @unpack bulk = model
    phiAB, phiCB, phiPB, phiMB = bulk
    x0[1] = phiAS/phiAB
    x0[2] = phiCS/phiCB
    x0[3] = phiPS/(1 - phiAS - phiCS)
    x0[4] = phiMS/(1 - phiAS - phiCS - phiPS)
    x0[5] = (phiAC - phiAB)/(1 - phiAB)
    x0[6] = (phiCC - phiCB)/(1 - phiAC - phiAB)
    x0[7] = phiPC/(1 - phiAC - phiCC)
    x0[8] = phiMC/(1 - phiAC - phiCC - phiPC)
    x0[9] = nu
    x0[10] = psi
    xs = [logscale(x0[1:9]), x0[10]]

end

function bndlunscale!(xs, model::TwoChainModel)

    x0 = [logunscale(xs[1:9]), x0[10]]
    phiAS = phiAB*x0[1]
    phiCS = phiCB*x0[2]
    phiPS = (1 - phiAS - phiCS)*x0[3]
    phiMS = (1 - phiAS - phiCS - phiPS)*x0[4]
    phiAC = phiAB + (1 - phiAB)*x0[5]
    phiCC = phiCB + (1 - phiAC - phiAB)*x0[6]
    phiPC = (1 - phiAC - phiCC)*x0[7]
    phiMS = (1 - phiAC - phiCC - phiPC)*x0[8]
    nu = x0[9]
    psi = x0[10]
    x = [phiAS, phiCS, phiPS, phiMS, phiAC, phiCC, phiPC, phiMC, nu, psi]

end

function bndlstate(x, model::TwoChainModel)
    @assert length(x) == 4 || length(x) == 10 "Invalid number of parameters for binodal state. Requires (4, 10) for $(typeof(model))."

    @unpack bulk = model
    phiAB, phiCB, phiPB, phiMB = bulk

    if length(x) == 4
        phiAC, phiCC, phiWC, nu = x
        phiWB = 1 - sum(bulk)
        @unpack omega, sig, zP, zM = model
        wA, wC, wP, wM = omega
        sigA, sigC = sig
        zwsum = zM*wP+zP*wM
        
        # Derive dense phase params
        ΔC = sigA*phiAC/wA - sigC*phiCC/zM
        if ΔC >= 0
            phiPC0 = ΔC*wP/zP
            phiMC0 = 0
        else
            phiPC0 = 0
            phiMC0 = -ΔC*wM/zM
        end
        phiPC = zM*wP/zwsum*(1 - (phiAC + phiCC + phiWC + phiPC0 + phiMC0)) + phiPC0
        phiMC = 1 - phiAC - phiCC - phiPC - phiWC
        
        # Derive supernatant phase params
        phiAS = (phiAB - nu*phiAC)/(1 - nu)
        phiCS = (phiCB - nu*phiCC)/(1 - nu)
        phiWS = (phiWB - nu*phiWC)/(1 - nu)
        ΔS = sigA*phiAS/wA - sigC*phiCS/zM
        if ΔS >= 0
            phiPS0 = ΔS*wP/zP
            phiMS0 = 0
        else
            phiPS0 = 0
            phiMs0 = -ΔS*wM/zM
        end
        phiPS = zM*wP/zwsum*(1 - (phiAS + phiCS + phiWS + phiPS0 + phiMS0)) + phiPS0
        phiMS = 1 - phiAS - phiCS - phiPS - phiWS

        # Update state struct
        bulk = @SVector [phiAB, phiCB, phiPB, phiMB]
        sup = @SVector [phiAS, phiCS, phiPS, phiMS]
        dense = @SVector [phiAC, phiCC, phiPC, phiMC]
        
        return BinodalState(bulk, sup, dense, nu)
    else
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
    @unpack bulk, zP, zM = model
    phiAB, phiCB, phiPB, phiMB = bulk
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
    F[3] = muS[3] - muC[3] - zP*psi/wP
    F[4] = muS[4] - muC[4] + zM*psi/wM
    F[5] = pS - pC
    F[6] = sigC*phiCS/wC - sigA*phiAS/wA + zP*phiPS/wP - zM*phiMS/wM
    F[7] = sigC*phiCC/wC - sigA*phiAC/wA + zP*phiPC/wP - zM*phiMC/wM
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
    z = model.z
    
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