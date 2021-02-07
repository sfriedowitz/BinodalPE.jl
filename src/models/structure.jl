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

struct ChainStructure{TC <: AbstractChainStructure}
	np :: Float64
	lp :: Float64
	b  :: Float64
	w  :: Float64
end

gchain(chain::ChainStructure{<:Union{PointLike,SmearedPoint}}, q) = 1.0

gchain(chain::ChainStructure{SphericalGlobule}, q) = gsphere(q * ((3/(4pi)) * chain.w * chain.np)^(1/3))

gchain(chain::ChainStructure{GaussianCoil}, q) = gcoil(q^2 * (chain.np * chain.b^2 / 6.0))

gchain(chain::ChainStructure{RodLike}, q) = gcoil(q * chain.np * chain.b)

gchain(chain::ChainStructure{<:Union{WormLike,AdaptiveChain}}, q) = gworm(q, chain.lp, chain.np, chain.b)