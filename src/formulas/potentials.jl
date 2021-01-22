#==============================================================================#
# Single polyion system
#==============================================================================#

muel_assoc(phi, vars, model::SinglePolyion{PointLike}) = model.lB*kbar(phi, vars, model)

muel_assoc(phi, vars, model::SinglePolyion{ExtendedPoint}) = model.lB*kbar(phi, vars, model)/(1 + kbar(phi, vars, model))

function muel_assoc(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    np = model.np
    b = model.b
    lB = model.lB

    # Derived parameters
    alpha = vars[1]
    sig = 1 - alpha

    # Integration function
    function integrand(q)
        gA = gchain(q, vars, model)
        smA, smP, smM = gamq.(q, model.smear)
        -(lB/pi)*(smP^2 + 2*np*sig*smA^2*gA)/(1 + ktilde2(q, phi, vars, model)/q^2)
    end

    pot, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER, maxevals = 100)
    return pot
end

#==============================================================================#
# Fully association coacervate
#==============================================================================#

function muel_assoc(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    wA, wC, wP, wM = model.omega
    nA, nC = model.np
    bA, bC = model.b
    lB = model.lB

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)

    # Integration functions
    function integrand(q)
        gA, gC = gchain(q, vars, model)
        smA, smC, smP, smM = gamq.(q, model.smear)

        num_ap = -(lB/pi)*(smP^2 + 2*nA*sigA*smA^2*gA)
        num_cm = -(lB/pi)*(smM^2 + 2*nC*sigC*smC^2*gC)
        num_bc  = -(2*lB/pi) * (nA*sigA*smA^2*gA + nC*sigC*smC^2*gC)
        numerators = @SVector [num_ap, num_cm, num_bc]

        numerators / (1 + ktilde2(q, phi, vars, model)/q^2)
    end

    pots, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER, maxevals = 100)
    return pots
end

function muel_assoc(phi, vars, model::AssociationCoacervate{PointLike})
    lB = model.lB
    kb = kbar(phi, vars, model)
    mu = lB*kb
    return @SVector [mu, mu, mu]
end

function muel_assoc(phi, vars, model::AssociationCoacervate{ExtendedPoint})
    lB = model.lB
    kb = kbar(phi, vars, model)
    mu = (lB*kb)/(1 + kb)
    return @SVector [mu, mu, mu]
end