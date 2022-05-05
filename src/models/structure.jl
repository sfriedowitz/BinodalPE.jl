 """
	abstract type AbstractChainStructure

An abstract representation of the polymer chain structure within a model.
A general random phase approximation expression is applied for the electrostatic free energy,
where a model's chain structure controls the polymer structure factor
used in the calculation.
"""
abstract type AbstractChainStructure end

#==============================================================================#

abstract type PointLike <: AbstractChainStructure end
abstract type ExtendedPoint <: AbstractChainStructure end
abstract type SmearedPoint <: AbstractChainStructure end
abstract type GaussianCoil <: AbstractChainStructure end
abstract type RodLike <: AbstractChainStructure end
abstract type WormLike <: AbstractChainStructure end
abstract type AdaptiveChain <: AbstractChainStructure end
abstract type SphericalGlobule <: AbstractChainStructure end

# temporary type for f_el using Edwards kernel and without smearing
abstract type EdwardsCoil <: AbstractChainStructure end
 
#==============================================================================#

struct ChainStructure{TC <: AbstractChainStructure, TF <: Real}
	w  :: Float64
	dp :: Float64
	b  :: Float64
	lp :: TF
end

gchain(q, chain::ChainStructure{<:Union{PointLike,SmearedPoint}}) = 1.0

gchain(q, chain::ChainStructure{SphericalGlobule}) = gsphere(q * ((3/(4pi)) * chain.w * chain.dp)^(1/3))

gchain(q, chain::ChainStructure{GaussianCoil}) = gcoil(q^2 * (chain.dp * chain.b^2 / 6.0))

gchain(q, chain::ChainStructure{RodLike}) = grod(q * chain.dp * chain.b)

gchain(q, chain::ChainStructure{<:Union{WormLike,AdaptiveChain}}) = gworm(q, chain.lp, chain.dp, chain.b)