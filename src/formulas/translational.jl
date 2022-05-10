#==============================================================================#
# Translational
#==============================================================================#

function ftranslational(phi, model::SinglePolyion)
    assoc = varsolve(phi, model)
    return ftranslational(phi, assoc, model)
end

function ftranslational(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    phiW = 1-sum(phi)
    alpha = vars[1]
    @unpack omega, dp, zP = model
    wA, wP, wM = omega
    nA = dp

    phiPF = phiP - alpha*phiA*wP/(wA*zP)

    return (phiA/wA/nA)*log(phiA) + (phiPF/wP)log(phiPF) + (phiM/wM)*log(phiM) + phiW*log(phiW)
end

function ftranslational(phi, model::SymmetricCoacervate)
    phiP, phiS = phi
    phiW = 1-sum(phi)
    @unpack omega, dp = model
    wP, wS = omega
    nP = dp
    return (phiP/wP/nP)*log(phiP) + (phiS/wS)*log(phiS) + phiW*log(phiW)
end

function ftranslational(phi, model::AsymmetricCoacervate)
    phiA, phiC, phiM, phiP = phi
    phiW = 1-sum(phi)
    @unpack omega, dp = model
    wA, wC, wP, wM = omega
    nA, nC = dp

    any(phi .< 0) && println(phi)

    ftrans = (phiA/wA/nA)*log(phiA) + (phiC/wC/nC)*log(phiC)
    ftrans += (phiP/wP)*log(phiP) + (phiM/wM)*log(phiM)
    ftrans += phiW*log(phiW)

    return ftrans
end

#==============================================================================#

function ftranslational(phi, model::AssociationCoacervate)
    assoc = varsolve(phi, model; model.vargs...)
    return ftranslational(phi, assoc, model)
end

function ftranslational(phi, assoc, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    phiW = 1-sum(phi)
    alphaAP, alphaCM, betaA, betaC = assoc
    @unpack omega, dp, zP, zM = model
    wA, wC, wP, wM = omega
    nA, nC = dp

    # Association adjusted fractions
    phiPF = phiP - (alphaAP*phiA*wP)/(wA*zP)
    phiMF = phiM - (alphaCM*phiC*wM)/(wC*zM)
    phiAB = phiA
    phiCB = phiC

    ftrans = (phiA/wA/nA)*log(phiA) + (phiC/wC/nC)*log(phiC)
    ftrans += (phiPF/wP)*log(phiPF) + (phiMF/wM)*log(phiMF)
    ftrans += phiW*log(phiW)

    return ftrans
end