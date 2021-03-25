#==============================================================================#
# Generic functional forms
#==============================================================================#

felectrostatic(phi, model::AbstractModel) = felectrostatic(phi, (), model)

function felectrostatic(phi, vars, model::AbstractModel)
    integrand(q) = (1/(4*π^2)) * q^2 * log(1 + kappa2(q, phi, vars, model)/q^2)
    sol, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return sol
end

function felectrostatic(phi, vars, model::AbstractModel{PointLike})
    kb = kbar(phi, vars, model)
    return -kb^3/(12π)
end

function felectrostatic(phi, vars, model::AbstractModel{ExtendedPoint})
    kb = kbar(phi, vars, model)
    return -1/(4π) * (log(1 + kb) - kb + (1/2)*kb^2)
end

#==============================================================================#
# Screening lengths
#==============================================================================#

kbar(phi, vars, model::AbstractModel) = notimpl("kbar", typeof(model))

function kbar(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    lB = model.lB

    # Derived parameters
    alpha = vars[1]    
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/wA

    return sqrt(4*π*lB*(phiPF/wP + phiM/wM + sigA*phiA/wA))
end

function kbar(phi, vars, model::SymmetricCoacervate)
    phiP, phiS = phi
    wP, wS = model.omega
    sig = model.sig
    lB = model.lB

    return sqrt(4*π*lB*(sig*phiP/wP + phiS/wS))
end

function kbar(phi, vars, model::AsymmetricCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    sigA, sigC = model.sig
    lB = model.lB

    return sqrt(4*π*lB*(sigA*phiA/wA + sigC*phiC/wC + phiP/wP + phiM/wM))
end

function kbar(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    nA, nC = model.dp
    lB = model.lB

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/wA
    phiMF = phiM - (alphaCM*phiC*wM)/wC

    return sqrt(4*π*lB*(sigA*phiA/wA + sigC*phiC/wC + phiPF/wP + phiMF/wM))
end

#==============================================================================#

kappa2(q, phi, vars, model::AbstractModel) = notimpl("kappa2", typeof(model))

function kappa2(q, phi, vars, model::SinglePolyion)
	phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    nA = model.dp
    lB = model.lB

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/wA

    smA, smP, smM = gamq.(q, model.smear)
    chain = chainstructs(model, vars)
    gA = gchain(q, chain)

	return (4*pi*lB)*(phiM*smM^2/wM + phiPF*smP^2/wP + phiA*nA*sigA^2*smA^2*gA/wA)
end

function kappa2(q, phi, vars, model::SymmetricCoacervate)
    phiP, phiS = phi
    wP, wS = model.omega
    sig = model.sig
    nP = model.dp
    lB = model.lB

    # Derived parameters
    chain = chainstructs(model, vars)
    gP = gchain(q, chain)
	smP, smS = gamq.(q, model.smear)

	return (4*pi*lB)*(phiS*smS^2/wS + phiP*nP*sig^2*smP^2*gP/wP)
end

function kappa2(q, phi, vars, model::AsymmetricCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    sigA, sigC = model.sig
    nA, nC = model.dp
    lB = model.lB

    # Derived parameters
    chains = chainstructs(model, vars)
    gA, gC = gchain.(q, chains)
    smA, smC, smP, smM = gamq.(q, model.smear)

	return (4*pi*lB)*(phiP*smP^2/wP + phiM*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end

function kappa2(q, phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    nA, nC = model.dp
    lB = model.lB

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/wA
    phiMF = phiM - (alphaCM*phiC*wM)/wC

    chains = chainstructs(model, vars)
    gA, gC = gchain.(q, chains)
    smA, smC, smP, smM = gamq.(q, model.smear)

	return (4*pi*lB)*(phiPF*smP^2/wP + phiMF*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end