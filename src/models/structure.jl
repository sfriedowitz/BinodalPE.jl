 """
	abstract type AbstractChainStructure

The representation of chain structure for use in defining the electrostatic free energy for a model.
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