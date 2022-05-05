#==============================================================================#
# Association combinatorics
#==============================================================================#

function fcombinatorial(phi, vars, model::SinglePolyion)
    phiA = phi[1]
    alpha = vars[1]
    @unpack omega = model
    wA = omega[1]

    return (phiA/wA)*(alpha*log(alpha) + (1-alpha)*log(1-alpha))
end

function fcombinatorial(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC = vars
    @unpack omega = model
    wA, wC, wP, wM = omega

    # Chain combinatorics terms
    fAP = (phiA/wA)*(alphaAP*log(alphaAP) + (1-alphaAP)*log(1-alphaAP))
    fCM = (phiC/wC)*(alphaCM*log(alphaCM) + (1-alphaCM)*log(1-alphaCM))

    fBA = ((1-alphaAP)*phiA/wA)*(betaA*log(betaA) + (1-betaA)*log(1-betaA))
    fBC = ((1-alphaCM)*phiC/wC)*(betaC*log(betaC) + (1-betaC)*log(1-betaC))

    # Bond volume factor
    fD = -(betaC*(1-alphaCM)*phiC/wC)*log(phiC*betaC*(1-alphaCM)*(wA+wC)/wC)

    return fAP + fCM + fBA + fBC + fD
end

#==============================================================================#
# Binding energy
#==============================================================================#

function fbinding(phi, vars, model::SinglePolyion)
    phiA = phi[1]
    alpha = vars[1]
    @unpack omega, dg = model
    wA = omega[1]

    return (phiA/wA)*alpha*dg
end

function fbinding(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC = vars
    @unpack dg, omega, zP, zM = model
    dgAP, dgCM, dgIP = dg
    wA, wC, wP, wM = omega

    return (alphaAP*phiA*dgAP)/(wA*zP) + (alphaCM*phiC*dgCM)/(wC*zM) + (betaC*(1-alphaCM)*phiC*dgIP)/wC
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