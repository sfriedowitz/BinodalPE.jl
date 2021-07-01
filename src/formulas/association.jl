#==============================================================================#
# Association combinatorics
#==============================================================================#

function fcombinatorial(phi, vars, model::SinglePolyion)
    phiA = phi[1]
    alpha = vars[1]
    wA = model.omega[1]

    return (phiA/wA)*(alpha*log(alpha) + (1-alpha)*log(1-alpha))
end

function fcombinatorial(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC = vars
    wA, wC, wP, wM = model.omega

    # Chain combinatorics terms
    fAP = (phiA/wA)*(alphaAP*log(alphaAP) + (1-alphaAP)*log(1-alphaAP))
    fCM = (phiC/wC)*(alphaCM*log(alphaCM) + (1-alphaCM)*log(1-alphaCM))

    fBA = ((1-alphaAP)*phiA/wA)*(betaA*log(betaA) + (1-betaA)*log(1-betaA))
    fBC = ((1-alphaCM)*phiC/wC)*(betaC*log(betaC) + (1-betaC)*log(1-betaC))

    # Bond volume factor
    fD = -(betaC*(1-alphaCM)*phiC/wC)*log(phiC*betaC*(1-alphaCM)*(wA+wC)/(wC*exp(1.0)))

    return fAP + fCM + fBA + fBC + fD
end

#==============================================================================#
# Binding energy
#==============================================================================#

function fbinding(phi, vars, model::SinglePolyion)
    phiA = phi[1]
    alpha = vars[1]
    wA = model.omega[1]
    dg = model.dg

    return (phiA/wA)*alpha*dg
end

function fbinding(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC = vars
    dgAP, dgCM, dgIP = model.dg
    wA, wC, wP, wM = model.omega

    return (alphaAP*phiA*dgAP)/(wA*z) + (alphaCM*phiC*dgCM)/(wC*z) + (betaC*(1-alphaCM)*phiC*dgIP)/wC
end

#==============================================================================#
# Self-complimentary binding
#==============================================================================#

function fselfcomp(phi, vars, model::SelfComplimentaryCoacervate)
    phiA = phi[1]
    tA = vars[end]
    wA = model.omega[1]
    fA = model.fA
    vAA = model.vAA
    dgAA = model.dgAA

    small = eps(eltype(phi))

    fchoice = (fA*phiA/wA) * (tA*log(tA + small) + (1-tA)*log(1-tA + small))
    fpair = -(fA*tA*phiA/2.0/wA) * log((fA*phiA*tA/exp(1.0))*(vAA/wA))
    fbind = (fA*tA*phiA/2.0/wA) * dgAA

    return fchoice + fpair + fbind
end
