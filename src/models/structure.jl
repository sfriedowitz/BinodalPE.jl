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
abstract type GaussianCoil <: AbstractChainStructure end
abstract type RodLike <: AbstractChainStructure end
abstract type WormLike <: AbstractChainStructure end
abstract type AdaptiveChain <: AbstractChainStructure end
abstract type SphericalGlobule <: AbstractChainStructure end