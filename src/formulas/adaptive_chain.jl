#==============================================================================#
# Adaptive chain structure free energy
#==============================================================================#

gamma2(lp, np, b = 1.0) = 2.0*(lp/b)*(1.0 - (lp/np/b)*(1.0 - exp(-np*b/lp)))

dgamma2(lp, np, b = 1.0) = (2.0/b)*(1.0 + exp(-np*b/lp)) + (4*lp/np/b^2)*(exp(-np*b/lp) - 1.0)

function dgworm(q, lp, np, b = 1.0)
    dcoil = -3*q*exp(-lp*q/2)*(6 + b*np*q*(2 + lp*q)) / (6 + b*lp*np*q^2)^2
    drod = q*exp(-lp*q/2)/(2*(1 + b*np*q/pi))
    return dcoil + drod
end

#==============================================================================#
# Single association polyanion
#==============================================================================#

function dfent_adaptive(vars, model::SinglePolyion{AdaptiveChain})
    lp = vars[2]
    np = model.np
    b = model.b

    g2 = gamma2(lp, np, b)
    dg2 = dgamma2(lp, np, b)
    dfent = -1.5*(dg2/g2 - dg2/(1 - g2/np))

    return dfent
end

function dfint_adaptive(phi, vars, model::SinglePolyion{AdaptiveChain})
    alpha, lp = vars
    phiA, phiP, phiM = phi
    wA, wP, wM = model.omega
    aA, aP, aM = model.smear
    np = model.np
    b = model.b
    lB = model.lB

    # Derived parameters
    sig = 1 - alpha
    phiPF = phiP - alpha*phiA*wP/wA

    ktilde2(q) = (4*pi*lB)*(phiPF*gamq(q,aP)^2/wP + phiM*gamq(q,aM)^2/wM + 
        phiA*np*sig^2*gamq(q,aA)^2*gworm(q,lp,np,b)/wA
    )
    spol(q) = sig^2 * np^2 * gamq(q,aA)^2 * gworm(q,lp,np,b)
    gscreen(q) = (4*pi*lB) / (q^2 + ktilde2(q))

    # Derivatives of the screened potentials and structure factors with lp
    dgscreen(q) = -(4*pi*lB)^2 * sig^2 * phiA*np/wA * gamq(q,aA)^2 * dgworm(q,lp,np,b) / (q^2 + ktilde2(q))^2
    dspol(q) = sig^2 * np^2 * gamq(q,aA)^2 * dgworm(q,lp,np,b)

    # From ZGW article on single-chain free energy minimization
    # dfint/dl = (1/4pi^2)∫dq q^2 d(G * S)/dl = (1/4pi^2)∫dq q^2[(dG/dL)*S + (dS/dL)*G]
    integrand(q) = (1/(4*pi^2)) * q^2 * (dgscreen(q) * spol(q) + gscreen(q) * dspol(q))

    dfint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return dfint
end

function dftot_adaptive(phi, vars, model::SinglePolyion{AdaptiveChain})
    dfent = dfent_adaptive(vars, model)
    dfint = dfint_adaptive(phi, vars, model)

    return dfent + dfint
end

#==============================================================================#
# Symmetric coacervate
#==============================================================================#

function fent_adaptive(vars, model::SymmetricCoacervate{AdaptiveChain})
    lp = vars[1]
    np = model.np
    b = model.b

    g2 = gamma2(lp, np, b)
    fent = -1.5*(np*log(1 - g2/np) + log(g2))

    return fent
end

function fint_adaptive(phi, vars, model::SymmetricCoacervate{AdaptiveChain})
    lp = vars[1]
    phiP, phiS = phi
    wP, wS = model.omega
    aP, aS = model.smear
    np = model.np
    sig = model.sig
    b = model.b
    lB = model.lB

    ktilde2(q) = (4*pi*lB)*(phiS*gamq(q,aS)^2/wS + phiP*np*sig^2*gamq(q,aP)^2*gworm(q,lp,np,b)/wP)
    spol(q) = sig^2 * np^2 * gamq(q,aP)^2 * gworm(q,lp,np,b)
    integrand(q) = (lB/pi) * q^2 * spol(q) / (q^2 + ktilde2(q))

    fint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return fint
end

function ftot_adaptive(phi, vars, model::SymmetricCoacervate{AdaptiveChain})
    fent = fent_adaptive(vars, model)
    fint = fint_adaptive(phi, vars, model)
    return fent + fint
end

#==============================================================================#

function dfent_adaptive(vars, model::SymmetricCoacervate{AdaptiveChain})
    lp = vars[1]
    np = model.np
    b = model.b

    g2 = gamma2(lp, np, b)
    dg2 = dgamma2(lp, np, b)
    dfent = -1.5*(dg2/g2 - dg2/(1 - g2/np))

    return dfent
end

function dfint_adaptive(phi, vars, model::SymmetricCoacervate{AdaptiveChain})
    lp = vars[1]
    phiP, phiS = phi
    wP, wS = model.omega
    aP, aS = model.smear
    np = model.np
    sig = model.sig
    b = model.b
    lB = model.lB

    ktilde2(q) = (4*pi*lB)*(phiS*gamq(q,aS)^2/wS + phiP*np*sig^2*gamq(q,aP)^2*gworm(q,lp,np,b)/wP)
    spol(q) = sig^2 * np^2 * gamq(q,aP)^2 * gworm(q,lp,np,b)
    gscreen(q) = (4*pi*lB) / (q^2 + ktilde2(q))

    # Derivatives of the screened potentials and structure factors with lp
    dgscreen(q) = -sig^2 * phiP*np/wP * gamq(q,aP)^2 * dgworm(q,lp,np,b) * gscreen(q)^2
    dspol(q) = sig^2 * np^2 * gamq(q,aP)^2 * dgworm(q,lp,np,b)

    # From ZGW article on single-chain free energy minimization
    # dfint/dl = (1/4pi^2)∫dq q^2 d(G * S)/dl = (1/4pi^2)∫dq q^2[(dG/dL)*S + (dS/dL)*G]
    integrand(q) = (1/(4*pi^2)) * q^2 * (dgscreen(q) * spol(q) + gscreen(q) * dspol(q))
    dfint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)

    return dfint
end

function dftot_adaptive(phi, vars, model::SymmetricCoacervate{AdaptiveChain})
    dfent = dfent_adaptive(vars, model)
    dfint = dfint_adaptive(phi, vars, model)
    return dfent + dfint
end

#==============================================================================#
# Asymmetric coacervate
#==============================================================================#

function dfent_adaptive(vars, model::AsymmetricCoacervate{AdaptiveChain})
    lpA, lpC = vars
    nA, nC = model.np
    bA, bC = model.b

    g2A = gamma2(lpA, nA, bA)
    g2C = gamma2(lpC, nC, bC)

    dg2A = dgamma2(lpA, nA, bA)
    dg2C = dgamma2(lpC, nC, bC)

    dfentA = -1.5*(dg2A/g2A - dg2A/(1 - g2A/nA))
    dfentC = -1.5*(dg2C/g2C - dg2C/(1 - g2C/nC))

    return (dfentA, dfentC)
end

function dfint_adaptive(phi, vars, model::AsymmetricCoacervate{AdaptiveChain})
    phiA, phiC, phiP, phiM = phi
    lpA, lpC = vars
    wA, wC, wP, wM = model.omega
    aA, aC, aP, aM = model.smear
    sigA, sigC = model.sig
    nA, nC = model.np
    bA, bC = model.b
    lB = model.lB

    ktilde2(q) = (4*pi*lB)*(phiP*gamq(q,aP)^2/wP + phiM*gamq(q,aM)^2/wM + 
        phiA*nA*sigA^2*gamq(q,aA)^2*gworm(q,lpA,nA,bA)/wA + 
        phiC*nC*sigC^2*gamq(q,aC)^2*gworm(q,lpC,nC,bC)/wC
    )
    spol(q) = @SVector [
        sigA^2 * nA^2 * gamq(q,aA)^2 * gworm(q,lpA,nA,bA), 
        sigC^2 * nC^2 * gamq(q,aC)^2 * gworm(q,lpC,nC,bC)
    ]
    gscreen(q) = (4*pi*lB) / (q^2 + ktilde2(q))

    # Derivatives of the screened potentials and structure factors with lp
    dgscreen(q) = @SVector [
        -(4*pi*lB)^2 * sigA^2 * phiA*nA/wA * gamq(q,aA)^2 * dgworm(q,lpA,nA,bA) / (q^2 + ktilde2(q))^2,
        -(4*pi*lB)^2 * sigC^2 * phiC*nC/wC * gamq(q,aC)^2 * dgworm(q,lpC,nC,bC) / (q^2 + ktilde2(q))^2
    ]
    dspol(q) = @SVector [
        sigA^2 * nA^2 * gamq(q,aA)^2 * dgworm(q,lpA,nA,bA),
        sigC^2 * nC^2 * gamq(q,aC)^2 * dgworm(q,lpC,nC,bC)
    ]

    # From ZGW article on single-chain free energy minimization
    # dfint/dl = (1/4pi^2)∫dq q^2 d(G * S)/dl = (1/4pi^2)∫dq q^2[(dG/dL)*S + (dS/dL)*G]
    integrand(q) = begin
        gsc = gscreen(q)
        dgsc = dgscreen(q)
        sp = spol(q)
        dsp = dspol(q)
        (1/(4*pi^2)) * q^2 * @SVector [
            dgsc[1] * sp[1] + gsc * dsp[1],
            dgsc[2] * sp[2] + gsc * dsp[2]
        ]
    end

    dfint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return dfint.data
end

function dftot_adaptive(phi, vars, model::AsymmetricCoacervate{AdaptiveChain})
    dfentA, dfentC = dfent_adaptive(vars, model)
    dfintA, dfintC = dfint_adaptive(phi, vars, model)

    return (dfentA + dfintA, dfentC + dfintC)
end

#==============================================================================#
# Association coacervate
#==============================================================================#

function fent_adaptive(vars, model::AssociationCoacervate{AdaptiveChain})
    _, _, _, _, lpA, lpC = vars
    nA, nC = model.np
    bA, bC = model.b

    g2A = gamma2(lpA, nA, bA)
    g2C = gamma2(lpC, nC, bC)

    fentA = -1.5*(nA*log(1 - g2A/nA) + log(g2A))
    fentC = -1.5*(nC*log(1 - g2C/nC) + log(g2C))

    return (fentA, fentC)
end

function fint_adaptive(phi, vars, model::AssociationCoacervate{AdaptiveChain})
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC, lpA, lpC = vars
    wA, wC, wP, wM = model.omega
    aA, aC, aP, aM = model.smear
    nA, nC = model.np
    bA, bC = model.b
    lB = model.lB

    # Derived parameters
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/wA
    phiMF = phiM - (alphaCM*phiC*wM)/wC

    ktilde2(q) = (4*pi*lB)*(phiPF*gamq(q,aP)^2/wP + phiMF*gamq(q,aM)^2/wM + 
        phiA*nA*sigA^2*gamq(q,aA)^2*gworm(q,lpA,nA,bA)/wA + 
        phiC*nC*sigC^2*gamq(q,aC)^2*gworm(q,lpC,nC,bC)/wC
    )
    spol(q) = @SVector [
        sigA^2 * nA^2 * gamq(q,aA)^2 * gworm(q,lpA,nA,bA), 
        sigC^2 * nC^2 * gamq(q,aC)^2 * gworm(q,lpC,nC,bC)
    ]
    integrand(q) = (lB/pi) * q^2 * spol(q) / (q^2 + ktilde2(q))

    fint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return fint.data
end

function ftot_adaptive(phi, vars, model::AssociationCoacervate{AdaptiveChain})
    fentA, fentC = fent_adaptive(vars, model)
    fintA, fintC = fint_adaptive(phi, vars, model)

    return (fentA + fintA, fentC + fintC)
end

#==============================================================================#

function dfent_adaptive(vars, model::AssociationCoacervate{AdaptiveChain})
    _, _, _, _, lpA, lpC = vars
    nA, nC = model.np
    bA, bC = model.b

    g2A = gamma2(lpA, nA, bA)
    g2C = gamma2(lpC, nC, bC)

    dg2A = dgamma2(lpA, nA, bA)
    dg2C = dgamma2(lpC, nC, bC)

    dfentA = -1.5*(dg2A/g2A - dg2A/(1 - g2A/nA))
    dfentC = -1.5*(dg2C/g2C - dg2C/(1 - g2C/nC))

    return (dfentA, dfentC)
end

function dfint_adaptive(phi, vars, model::AssociationCoacervate{AdaptiveChain})
    phiA, phiC, phiP, phiM = phi
    alphaAP, alphaCM, betaA, betaC, lpA, lpC = vars
    wA, wC, wP, wM = model.omega
    aA, aC, aP, aM = model.smear
    nA, nC = model.np
    bA, bC = model.b
    lB = model.lB

    # Derived parameters
    sigA = (1-alphaAP)*(1-betaA)
    sigC = (1-alphaCM)*(1-betaC)
    phiPF = phiP - (alphaAP*phiA*wP)/wA
    phiMF = phiM - (alphaCM*phiC*wM)/wC

    ktilde2(q) = (4*pi*lB)*(phiPF*gamq(q,aP)^2/wP + phiMF*gamq(q,aM)^2/wM + 
        phiA*nA*sigA^2*gamq(q,aA)^2*gworm(q,lpA,nA,bA)/wA + 
        phiC*nC*sigC^2*gamq(q,aC)^2*gworm(q,lpC,nC,bC)/wC
    )
    spol(q) = @SVector [
        sigA^2 * nA^2 * gamq(q,aA)^2 * gworm(q,lpA,nA,bA), 
        sigC^2 * nC^2 * gamq(q,aC)^2 * gworm(q,lpC,nC,bC)
    ]
    gscreen(q) = (4*pi*lB) / (q^2 + ktilde2(q))

    # Derivatives of the screened potentials and structure factors with lp
    dgscreen(q) = @SVector [
        -(4*pi*lB)^2 * sigA^2 * phiA*nA/wA * gamq(q,aA)^2 * dgworm(q,lpA,nA,bA) / (q^2 + ktilde2(q))^2,
        -(4*pi*lB)^2 * sigC^2 * phiC*nC/wC * gamq(q,aC)^2 * dgworm(q,lpC,nC,bC) / (q^2 + ktilde2(q))^2
    ]
    dspol(q) = @SVector [
        sigA^2 * nA^2 * gamq(q,aA)^2 * dgworm(q,lpA,nA,bA),
        sigC^2 * nC^2 * gamq(q,aC)^2 * dgworm(q,lpC,nC,bC)
    ]

    # From ZGW article on single-chain free energy minimization
    # dfint/dl = (1/4pi^2)∫dq q^2 d(G * S)/dl = (1/4pi^2)∫dq q^2[(dG/dL)*S + (dS/dL)*G]
    integrand(q) = begin
        gsc = gscreen(q)
        dgsc = dgscreen(q)
        sp = spol(q)
        dsp = dspol(q)
        (1/(4*pi^2)) * q^2 * @SVector [
            dgsc[1] * sp[1] + gsc * dsp[1],
            dgsc[2] * sp[2] + gsc * dsp[2]
        ]
    end

    dfint, _ = quadgk(integrand, 0.0, Inf, order = QGK_ORDER)
    return dfint.data
end

function dftot_adaptive(phi, vars, model::AssociationCoacervate{AdaptiveChain})
    dfentA, dfentC = dfent_adaptive(vars, model)
    dfintA, dfintC = dfint_adaptive(phi, vars, model)

    return (dfentA + dfintA, dfentC + dfintC)
end