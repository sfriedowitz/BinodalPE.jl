#==============================================================================#
# Chi-Interaction
#==============================================================================#

function fchi(phi, model::Union{SinglePolyion,SymmetricCoacervate})
	phiP = phi[1]
	phiW = 1-sum(phi)
	return model.chi * phiP * phiW
end

function fchi(phi, model::AsymmetricCoacervateModel)
	phiA, phiC = phi[1], phi[2]
	phiW = 1-sum(phi)
	chiAW, chiCW, chiAC = model.chi

	return chiAW*phiA*phiW + chiCW*phiC*phiW + chiAC*phiA*phiC
end