"""
    AssociationCoacervate{TC <: AbstractChainStucture} <: AbstractModel{TC}

A coacervate model with treament of reversible binding 
for strongly dissociating polyelectrolytes.

At a given composition, the extent of each reaction is determined with the `association_solve` method.
These reversible reactions represent:

1) Cation binding on the polyanions
2) Anion binding on the polycations
3) Formation of inter-chain ion-pairs between the polyions

References:
* Salehi, A.; Larson, R. G. Macromolecules 2016, 49 (24), 9706–9719.
* Friedowitz, S.; Salehi, A.; Larson, R. G.; Qin, J. J. Chem. Phys. 2018, 149 (16), 163335.
"""
mutable struct AssociationCoacervate{TC <: AbstractChainStructure} <: AbstractModel{TC}
    bulk     :: MVector{4,Float64}
    omega    :: SVector{4,Float64}
    smear    :: SVector{4,Float64}
    sig      :: SVector{2,Float64}
    dg       :: SVector{3,Float64}
    chi      :: SVector{3,Float64}
    np       :: SVector{2,Float64}
    b        :: SVector{2,Float64}
    lp       :: SVector{2,Float64}
    lB       :: Float64
    vargs    :: Dict{Symbol,Any}
end

function AssociationCoacervate(; structure::Type{<:AbstractChainStructure}, kwargs...)
    omega = get(kwargs, :omega, [1.0, 1.0, 1.0, 1.0])
    sig = get(kwargs, :sig, [1.0, 1.0])
    dg = get(kwargs, :dg, [0.0, 0.0, 0.0])
    chi = get(kwargs, :chi, [0.0, 0.0, 0.0])
    np = get(kwargs, :np, [100.0, 100.0])
    b = get(kwargs, :b, [1.0, 1.0])
    lp = get(kwargs, :lp, [1.0, 1.0])
    lB = get(kwargs, :lB, lBbar)
    vargs = get(kwargs, :vargs, Dict())

    smear = 0.5 .* omega .^ (1/3)
    model = AssociationCoacervate{structure}(zeros(4), omega, smear, sig, dg, chi, np, b, lp, lB, vargs)
    return model
end

#==============================================================================#

function ftotal(phi, model::AssociationCoacervate{TC}) where TC
    # Solve for association fractions
    assoc = zeros(eltype(phi), 4)
    if TC == AdaptiveChain
        push!(assoc, model.b[1]/2, model.b[2]/2)
    end
    try
        assoc :: Vector{eltype(phi)} = varsolve(phi, model; model.vargs...)
    catch e
        @warn "Association solver failure." maxlog = 1
        nothing
    end

    ftrans = ftranslational(phi, assoc, model)
    fcomb = fcombinatorial(phi, assoc, model)
    finf = fbinding(phi, assoc, model)
    ffh = fchi(phi, model)
    fel = felectrostatic(phi, assoc, model)

    return ftrans + fcomb + finf + ffh + fel
end

#==============================================================================#
# Variational solver functions
#==============================================================================#

varscale(x, model::AssociationCoacervate) = [log(v/(1-v)) for v in x]

varunscale(x, model::AssociationCoacervate) = [exp(v)/(1 + exp(v)) for v in x]

function varscale(x, model::AssociationCoacervate{AdaptiveChain})
    new = similar(x)
    for i = 1:4; new[i] = log(x[i]/(1-x[i])); end
    new[5] = log(x[4])
    new[6] = log(x[5])
    return new
end

function varunscale(x, model::AssociationCoacervate{AdaptiveChain})
    new = similar(x)
    for i = 1:4; new[i] = exp(x[i])/(1+exp(x[i]));  end
    new[5] = exp(x[4])
    new[6] = exp(x[5])
    return new
end

other_beta(phi, omega, vars) = other_beta(phi[1], phi[2], omega[1], omega[2], vars[1], vars[2], vars[3])

other_beta(phiA, phiC, wA, wC, alphaAP, alphaCM, betaA) = phiA*betaA*(1-alphaAP)/wA * (wC/(phiC*(1-alphaCM)))

function valid_variational(xs, phi, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    
    alphaAP = exp(xs[1]) / (1 + exp(xs[1]))
    alphaCM = exp(xs[2]) / (1 + exp(xs[2]))
    
    phiPF = phiP - alphaAP*phiA*wP/wA
    phiMF = phiM - alphaCM*phiC*wM/wC
    
    return (0.0 < phiPF < 1.0) && (0.0 < phiMF < 1.0)   
end

function varinit(phi, model::AssociationCoacervate{TC}) where TC
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    dgAP, dgCM, dgIP = model.dg
    bA, bC = model.b
    nA, nC = model.np

    mu_ap, mu_cm, mu_ip = muel_association(phi, [0.5, 0.5, 0.5, 0.5, bA, bC], model)
    kAP = exp(-dgAP - mu_ap + 1)
    kCM = exp(-dgCM - mu_cm + 1)
    kIP = exp(-dgIP - mu_ip + 1)

    # Guess for alphaAP
    if phiA < 1e-10
        alphaAP0 = kAP*phiP/(1 + kAP*phiP)
    else
        innerAP = (wA + kAP*phiP*wA)^2 - 2*kAP*phiA*(kAP*phiP - 1)*wA*wP + (kAP*phiA*wP)^2
        quadAP = innerAP < 0 ? 1 : sqrt(innerAP)
        alphaAP0 = wA/(2*kAP*phiA*wP) + (phiP*wA + phiA*wP)/(2*phiA*wP) - quadAP/(2*kAP*phiA*wP)

        if !(0.0 < alphaAP0 < 1.0)
            alphaAP0 = wA/(2*kAP*phiA*wP) + (phiP*wA + phiA*wP)/(2*phiA*wP) + quadAP/(2*kAP*phiA*wP)
        end
    end

    # Guess for alphaCM
    if phiC < 1e-10
        alphaCM0 = kCM*phiM/(1 + kCM*phiM)
    else
        innerCM = (wC + kCM*phiM*wC)^2 - 2*kCM*phiC*(kCM*phiM - 1)*wC*wM + (kCM*phiC*wM)^2
        quadCM = innerCM < 0 ? 1 : sqrt(innerCM)
        alphaCM0 = wC/(2*kCM*phiC*wM) + (phiM*wC + phiC*wM)/(2*phiC*wM) - quadCM/(2*kCM*phiC*wM)
        if !(0.0 < alphaCM0 < 1.0)
            alphaCM0 = wC/(2*kCM*phiC*wM) + (phiM*wC + phiC*wM)/(2*phiC*wM) + quadCM/(2*kCM*phiC*wM)
        end
    end

    alphaAP0 = clamp(alphaAP0, 1e-14, 0.9999)
    alphaCM0 = clamp(alphaCM0, 1e-14, 0.9999)

    # Ion pairing extents
    betaA = (1-alphaAP0)*phiA/wA
    betaC = (1-alphaCM0)*phiC/wC
    
    M = min(betaA, betaC) / max(betaA, betaC)
    G = (wA+wC)*max(betaA, betaC)
    if M < 1e-4
        beta = kIP*G / (1 + kIP*G)
    elseif kIP*G < 1e-6
        beta = sqrt(kIP*G)
    else
        beta = (1/(2*M))*(1 + M + (1/kIP/G) + sqrt((1 + M + 1/kIP/G)^2 - 4*M))
    end
    if !(0.0 < beta < 1.0)
        beta = (1/(2*M))*(1 + M + (1/kIP/G) - sqrt((1 + M + 1/kIP/G)^2 - 4*M))
    end
    
    if betaC <= betaA
        betaA0 = beta * betaC / betaA
        betaC0 = beta
    else
        betaA0 = beta
        betaC0 = beta * betaA / betaC
    end

    init = [alphaAP0, alphaCM0, betaA0, betaC0]
    for i in eachindex(init)
        if init[i] <= 0.0 || init[i] >= 1.0 || isnan(init[i])
            init[i] = 0.5
        end
    end

    if TC == AdaptiveChain
        push!(init, bA/2, bC/2)
    end

    return init
end

#==============================================================================#

function varf!(F, x, phi, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    dgAP, dgCM, dgIP = model.dg
   
    vars = varunscale(x, model)
    alphaAP, alphaCM, betaA, betaC = vars
    #betaC = other_beta(phiA, phiC, wA, wC, alphaAP, alphaCM, betaA)

    phiPF = phiP - alphaAP*phiA*wP/wA
    phiMF = phiM - alphaCM*phiC*wM/wC
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    
    # Electrostatic exchange potentials
    mu_ap, mu_cm, mu_bc = muel_association(phi, vars, model)
    
    F[1] = log(alphaAP/(sigA*phiPF)) + dgAP + mu_ap - 1
    F[2] = log(alphaCM/(sigC*phiMF)) + dgCM + mu_cm - 1
    F[3] = log(exp(1.0) * (betaC/(1-betaC)) * (wA/(wA + wC)) / (sigA*phiA)) + dgIP + mu_bc - 1
    F[4] = phiA*betaA*(1-alphaAP)/wA - phiC*betaC*(1-alphaCM)/wC

    return nothing
end

function varf!(F, x, phi, model::AssociationCoacervate{AdaptiveChain})
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    dgAP, dgCM, dgIP = model.dg

    vars = varunscale(x, model)
    alphaAP, alphaCM, betaA, betaC, lpA, lpC = vars
    #betaC = other_beta(phiA, phiC, wA, wC, alphaAP, alphaCM, betaA)

    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - alphaAP*phiA*wP/wA
    phiMF = phiM - alphaCM*phiC*wM/wC
    
    # Electrostatic exchange potentials
    mu_ap, mu_cm, mu_bc = muel_association(phi, vars, model)
    dfA, dfC = dftot_adaptive(phi, vars, model)
    
    F[1] = log(alphaAP/(sigA*phiPF)) + dgAP + mu_ap - 1
    F[2] = log(alphaCM/(sigC*phiMF)) + dgCM + mu_cm - 1
    F[3] = log((betaC*exp(1.0)*wA)/(1-betaC)/(sigA*phiA)/(wA+wC)) + dgIP + mu_bc - 1
    F[4] = phiA*betaA*(1-alphaAP)/wA - phiC*betaC*(1-alphaCM)/wC
    F[5] = dfA
    F[6] = dfC

    return nothing
end