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
 
#==============================================================================#

struct ChainStructure{TC <: AbstractChainStructure, TF <: Real}
	w  :: Float64
	dp :: Float64
	b  :: Float64
	lp :: TF
end

gchain(chain::ChainStructure{<:Union{PointLike,SmearedPoint}}, q) = 1.0

gchain(chain::ChainStructure{SphericalGlobule}, q) = gsphere(q * ((3/(4pi)) * chain.w * chain.dp)^(1/3))

gchain(chain::ChainStructure{GaussianCoil}, q) = gcoil(q^2 * (chain.dp * chain.b^2 / 6.0))

gchain(chain::ChainStructure{RodLike}, q) = grod(q * chain.dp * chain.b)

gchain(chain::ChainStructure{<:Union{WormLike,AdaptiveChain}}, q) = gworm(q, chain.lp, chain.dp, chain.b)