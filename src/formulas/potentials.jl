 #==============================================================================#
# Single polyion system
#==============================================================================#

muel_association(phi, vars, model::SinglePolyion{PointLike}) = model.lB*kbar(phi, vars, model)

muel_association(phi, vars, model::SinglePolyion{ExtendedPoint}) = model.lB*kbar(phi, vars, model)/(1 + kbar(phi, vars, model))

function muel_association(phi, vars, model::SinglePolyion)
    phiA, phiP, phiM = phi
    @unpack omega, dp, lB, zP = model
    wA, wP, wM = omega
    nA = dp

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    chain = chainstructs(model, vars)

    # Integration function
    function integrand(q)
        gA = gchain(q, chain)
        smA, smP, smM = gamq.(q, model.smear)
        -(lB/pi)*(zP*smP^2 + 2*nA*sigA*smA^2*gA)/(1 + kappa2(q, phi, vars, model)/q^2)
    end

    # # Alternate integration function to possibly avoid precision issues at low q
    # function integrand(q)
    #     gA = gchain(q, chain)
    #     smA, smP, smM = gamq.(q, model.smear)
    #     q2 = q*q
    #     int = -q2*(lB/pi)*(z*smP^2 + 2*nA*sigA*smA^2*gA)/(q2 + kappa2(q, phi, vars, model))
    #     # diagnostic output
    #     isnan(int) && println("NaN at alpha = $alpha, q = $q, phi = $phi")
    #     return int
    # end    

    pot, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER, maxevals = 100)
    return pot
end

function muel_association(phi, vars, model::SinglePolyion{SmearedPoint})
    phiA, phiP, phiM = phi
    @unpack omega, dp, lB, zP = model
    wA, wP, wM = omega
    nA = dp

    # Derived parameters
    alpha = vars[1]
    sigA = 1 - alpha
    chain = chainstructs(model, vars)

    # Integration function
    function integrand(q)
        smA, smP, smM = gamq.(q, model.smear)
        -(lB/pi)*(zP*smP^2 + smA^2)/(1 + kappa2(q, phi, vars, model)/q^2)
    end

    pot, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER, maxevals = 100)
    return pot
end

#==============================================================================#
# Fully association coacervate
#==============================================================================#

function muel_association(phi, vars, model::AssociationCoacervate)
    phiA, phiC, phiP, phiM = phi
    @unpack omega, dp, lB, zP, zM = model
    wA, wC, wP, wM = omega
    nA, nC = dp

    # Derived parameters
    alphaAP, alphaCM, betaA, betaC = vars
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)

    chains = chainstructs(model, vars)

    # Integration functions
    function integrand(q)
        gA, gC = gchain.(q, chains)
        smA, smC, smP, smM = gamq.(q, model.smear)

        num_ap = -(lB/pi)*(zP*smP^2 + 2*nA*sigA*smA^2*gA)
        num_cm = -(lB/pi)*(zM*smM^2 + 2*nC*sigC*smC^2*gC)
        num_bc  = -(2*lB/pi) * (nA*sigA*smA^2*gA + nC*sigC*smC^2*gC)
        numerators = @SVector [num_ap, num_cm, num_bc]

        numerators / (1 + kappa2(q, phi, vars, model)/q^2)
    end

    pots, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER, maxevals = 100)
    return pots
end

function muel_association(phi, vars, model::AssociationCoacervate{PointLike})
    lB = model.lB
    kb = kbar(phi, vars, model)
    mu = lB*kb
    return @SVector [mu, mu, mu]
end

function muel_association(phi, vars, model::AssociationCoacervate{ExtendedPoint})
    lB = model.lB
    kb = kbar(phi, vars, model)
    mu = (lB*kb)/(1 + kb)
    return @SVector [mu, mu, mu]
end