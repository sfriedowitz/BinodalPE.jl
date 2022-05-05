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


# temporary type for f_el using Edwards kernel and without smearing
function felectrostatic(phi, vars, model::SymmetricCoacervate{EdwardsCoil})
    @unpack omega, sig, dp, lB, zP, zM, b = model
    wP, wS = omega
    phiP, phiS = phi
    phiP /= wP
    phiS /= wS
    sig = sig

    # Derived parameters
    p = (3*pi*lB*sig^2/b^2)^(3/4)/(3*pi)
    q = sqrt(pi*lB*b^2/(3*sig^2))
    s = zP*zM*q*phiS/sqrt(phiP)


    return 2*p*phiP^(3/4)*(1 - s)*sqrt(2 + s)
end

#==============================================================================#
# Screening lengths
#==============================================================================#

kbar(phi, vars, model::AbstractModel) = notimpl("kbar", typeof(model))

function kbar(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    @unpack omega, lB, zP, zM = model
    wA, wP, wM = omega

    # Derived parameters
    alpha = vars[1]    
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/(wA*zP)

    return sqrt(4*π*lB*(zP^2*phiPF/wP + zM^2*phiM/wM + sigA*phiA/wA))
end

function kbar(phi, vars, model::SymmetricCoacervate)
    phiP, phiS = phi
    @unpack omega, sig, lB, zP, zM = model
    wP, wS = omega

    return sqrt(4*π*lB*(sig*phiP/wP + zP*zM*phiS/wS))
end

function kbar(phi, vars, model::AsymmetricCoacervate)
    phiA, phiC, phiP, phiM = phi
    @unpack omega, sig, lB, zP, zM = model
    wA, wC, wP, wM = omega
    sigA, sigC = sig

    return sqrt(4*π*lB*(sigA*phiA/wA + sigC*phiC/wC + zP^2*phiP/wP + zM^2*phiM/wM))
end

function kbar(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    @unpack omega, dp, lB, zP, zM = model
    wA, wC, wP, wM = omega
    nA, nC = dp

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/(wA*zP)
    phiMF = phiM - (alphaCM*phiC*wM)/(wC*zM)

    return sqrt(4*π*lB*(sigA*phiA/wA + sigC*phiC/wC + zP^2*phiPF/wP + zM^2*phiMF/wM))
end

#==============================================================================#

kappa2(q, phi, vars, model::AbstractModel) = notimpl("kappa2", typeof(model))

function kappa2(q, phi, vars, model::SinglePolyion)
	phiA, phiP, phiM = phi
    @unpack omega, dp, lB, zP, zM, smear = model
    wA, wP, wM = omega
    nA = dp

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/(wA*zP)

    smA, smP, smM = gamq.(q, smear)
    chain = chainstructs(model, vars)
    gA = gchain(q, chain)

	return (4*pi*lB)*(zM^2*phiM*smM^2/wM + zP^2*phiPF*smP^2/wP + phiA*nA*sigA^2*smA^2*gA/wA)
end

function kappa2(q, phi, vars, model::SinglePolyion{SmearedPoint})
    phiA, phiP, phiM = phi
    @unpack omega, dp, lB, zP, zM, smear = model
    wA, wP, wM = omega
    nA = dp

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/(wA*zP)

    smA, smP, smM = gamq.(q, smear)

    return (4*pi*lB)*(zM^2*phiM*smM^2/wM + zP^2*phiPF*smP^2/wP + phiA*sigA*smA^2/wA)
end

function kappa2(q, phi, vars, model::SymmetricCoacervate)
    phiP, phiS = phi
    @unpack omega, sig, dp, lB, zP, zM, smear = model
    wP, wS = omega
    sig = sig
    nP = dp

    # Derived parameters
    chain = chainstructs(model, vars)
    gP = gchain(q, chain)
	smP, smS = gamq.(q, smear)

	return (4*pi*lB)*(zP*zM*phiS*smS^2/wS + phiP*nP*sig^2*smP^2*gP/wP)
end

function kappa2(q, phi, vars, model::SymmetricCoacervate{SmearedPoint})
    phiP, phiS = phi
    @unpack omega, sig, dp, lB, zP, zM, smear = model
    wP, wS = omega
    sig = sig
    nP = dp

    # Derived parameters
    smP, smS = gamq.(q, smear)

    return (4*pi*lB)*(zP*zM*phiS*smS^2/wS + phiP*sig*smP^2/wP)
end

function kappa2(q, phi, vars, model::AsymmetricCoacervate)
    phiA, phiC, phiP, phiM = phi
    @unpack omega, sig, dp, lB, zP, zM, smear = model
    wA, wC, wP, wM = omega
    sigA, sigC = sig
    nA, nC = dp

    # Derived parameters
    chains = chainstructs(model, vars)
    gA, gC = gchain.(q, chains)
    smA, smC, smP, smM = gamq.(q, smear)

	return (4*pi*lB)*(zP^2*phiP*smP^2/wP + zM^2*phiM*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end

function kappa2(q, phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    @unpack omega, dp, lB, zP, zM, smear = model
    wA, wC, wP, wM = omega
    nA, nC = dp

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/(wA*zP)
    phiMF = phiM - (alphaCM*phiC*wM)/(wC*zM)

    chains = chainstructs(model, vars)
    gA, gC = gchain.(q, chains)
    smA, smC, smP, smM = gamq.(q, smear)

	return (4*pi*lB)*(zP^2*phiPF*smP^2/wP + zM^2*phiMF*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end