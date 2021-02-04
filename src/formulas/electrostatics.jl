#==============================================================================#
# Generic function for all
#==============================================================================#

felectrostatic(phi, model::AbstractModel) = felectrostatic(phi, (), model)

function felectrostatic(phi, vars, model::AbstractModel)
    integrand(q) = (1/(4*π^2)) * q^2 * log(1 + ktilde2(q, phi, vars, model)/q^2)
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
# Chain structure form factors
#==============================================================================#

# Single polyion and symmetric coacervate

const OneChainModel{TC} = Union{SinglePolyion{TC}, SymmetricCoacervate{TC}}
gchain(q, vars, model::OneChainModel{SmearedPoint}) = 1.0
gchain(q, vars, model::OneChainModel{SphericalGlobule}) = gsphere(q*((3/(4pi))*model.omega[1]*model.np)^(1/3))
gchain(q, vars, model::OneChainModel{GaussianCoil}) = gcoil(q^2*(model.np*model.b^2/6.0)) 
gchain(q, vars, model::OneChainModel{RodLike}) = grod(q*model.np*model.b) 
gchain(q, vars, model::OneChainModel{WormLike}) = gworm(q, model.lp, model.np, model.b)
gchain(q, vars, model::OneChainModel{AdaptiveChain}) = gworm(q, vars[end], model.np, model.b)

# Asymmetric coacervate models

gchain(q, vars, model::AsymmetricCoacervateModel{SmearedPoint}) = (1.0, 1.0)

function gchain(q, vars, model::AsymmetricCoacervateModel{SphericalGlobule})
    wA, wC = model.omega[1], model.omega[2]
    nA, nC = model.np
    rA = ((3/(4pi))*nA*wA)^(1/3)
    rC = ((3/(4pi))*nC*wC)^(1/3)
    return (gsphere(q*rA), gsphere(q*rC))
end

function gchain(q, vars, model::AsymmetricCoacervateModel{GaussianCoil})
    nA, nC = model.np
    bA, bC = model.b
    return (gcoil(q^2*nA*bA^2/6.0), gcoil(q^2*nC*bC^2/6.0))
end

function gchain(q, vars, model::AsymmetricCoacervateModel{RodLike})
    nA, nC = model.np
    bA, bC = model.b
    return (grod(q*nA*bA), grod(q*nC*bC))
end

function gchain(q, vars, model::AsymmetricCoacervateModel{WormLike})
    lpA, lpC = model.lp
    nA, nC = model.np
    bA, bC = model.b
    return (gworm(q,lpA,nA,bA), gworm(q,lpC,nC,bC))
end

function gchain(q, vars, model::AsymmetricCoacervateModel{AdaptiveChain})
    lpA, lpC = vars[end-1], vars[end]
    nA, nC = model.np
    bA, bC = model.b
    return (gworm(q,lpA,nA,bA), gworm(q,lpC,nC,bC))
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
    nA, nC = model.np
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

ktilde2(q, phi, vars, model::AbstractModel) = notimpl("ktilde2", typeof(model))

function ktilde2(q, phi, vars, model::SinglePolyion)
	phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    nA = model.np
    lB = model.lB

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/wA
	
	gA = gchain(q, vars, model)
	smA, smP, smM = gamq.(q, model.smear)

	return (4*pi*lB)*(phiM*smM^2/wM + phiPF*smP^2/wP + phiA*nA*sigA^2*smA^2*gA/wA)
end

function ktilde2(q, phi, vars, model::SymmetricCoacervate)
    phiP, phiS = phi
    wP, wS = model.omega
    sig = model.sig
    np = model.np
    lB = model.lB

    # Derived parameters
    gP = gchain(q, vars, model)
	smP, smS = gamq.(q, model.smear)

	return (4*pi*lB)*(phiS*smS^2/wS + phiP*np*sig^2*smP^2*gP/wP)
end

function ktilde2(q, phi, vars, model::AsymmetricCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    sigA, sigC = model.sig
    nA, nC = model.np
    lB = model.lB

    # Derived parameters
    gA, gC = gchain(q, vars, model)
	smA, smC, smP, smM = gamq.(q, model.smear)

	return (4*pi*lB)*(phiP*smP^2/wP + phiM*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end

function ktilde2(q, phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    nA, nC = model.np
    lB = model.lB

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/wA
    phiMF = phiM - (alphaCM*phiC*wM)/wC

    gA, gC = gchain(q, vars, model)
	smA, smC, smP, smM = gamq.(q, model.smear)

	return (4*pi*lB)*(phiPF*smP^2/wP + phiMF*smM^2/wM + phiA*nA*sigA^2*smA^2*gA/wA + phiC*nC*sigC^2*smC^2*gC/wC)
end