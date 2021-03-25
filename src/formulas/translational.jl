#==============================================================================#
# Translational
#==============================================================================#

function ftranslational(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    phiW = 1-sum(phi)
    alpha = vars[1]
    wA, wP, wM = model.omega
    nA = model.dp

    phiAB = phiA + alpha*phiA*wP/wA
    phiPF = phiP - alpha*phiA*wP/wA

    return (phiAB/wA/nA)*log(phiAB) + (phiPF/wP)log(phiPF) + (phiM/wM)*log(phiM) + phiW*log(phiW)
end

function ftranslational(phi, model::SymmetricCoacervate)
    phiP, phiS = phi
    phiW = 1-sum(phi)
    wP, wS = model.omega
    nP = model.dp
    return (phiP/wP/nP)*log(phiP) + (phiS/wS)*log(phiS) + phiW*log(phiW)
end

function ftranslational(phi, model::AsymmetricCoacervate)
    phiA, phiC, phiM, phiP = phi
    phiW = 1-sum(phi)
    wA, wC, wP, wM = model.omega
    nA, nC = model.dp

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
    wA, wC, wP, wM = model.omega
    nA, nC = model.dp

    # Association adjusted fractions
    phiPF = phiP - (alphaAP*phiA*wP)/wA    
    phiMF = phiM - (alphaCM*phiC*wM)/wC
    phiAB = phiA + (alphaAP*phiA*wP)/wA
    phiCB = phiC + (alphaCM*phiC*wM)/wC

    ftrans = (phiAB/wA/nA)*log(phiAB) + (phiCB/wC/nC)*log(phiCB)
    ftrans += (phiPF/wP)*log(phiPF) + (phiMF/wM)*log(phiMF)
    ftrans += phiW*log(phiW)

    return ftrans
end