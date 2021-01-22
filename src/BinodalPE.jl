module BinodalPE

#==============================================================================#
# Imports
#==============================================================================#

using Printf
using JSON
using RecipesBase
using DelimitedFiles

using QuadGK
using StaticArrays
using Interpolations
using LinearAlgebra

using ForwardDiff
import ForwardDiff: Chunk, HessianConfig, JacobianConfig

using NLsolve
using Optim

import Distances: sqeuclidean
import SpecialFunctions: sinint

#==============================================================================#
# Exported types
#==============================================================================#

# Types
export AbstractChainStructure, PointLike, ExtendedPoint, GaussianCoil, RodLike, WormLike, AdaptiveChain, SphericalGlobule
export AbstractModel, AsymmetricCoacervateModel
export SinglePolyion, SymmetricCoacervate, AsymmetricCoacervate, AssociationCoacervate, SelfComplimentaryCoacervate

# Solve binodals
export BinodalState, BinodalData, BinodalResults
export bndlminx, bndlsolvex, bndlstate, bndlscale, bndlunscale
export bndlf, bndlf!, bndlj, bndlj!, bndlfj!, bndlsolve
export bndlg, bndlproblem, bndlminimize
export savebndl, readbndl, add_state!

# Model methods
export swap, swap!, valid, set_bulk, newstate
export free_energy, fideal, fexcess, f2total, f3total, mutotal, muexcess, pressure, mupressure
export varinit, varsolve, varscale, varunscale, varf, varf!, varj!, varfj!

# Some utilities
export gcoil, grod, gworm, gchain
export tophi, toconc, neutralbulk, differencebulk, asypolyion

#==============================================================================#
# Constants
#==============================================================================#

const EPSF64 = eps(Float64)
const QGK_ORDER = 50

#==============================================================================#
# Included files
#==============================================================================#

include("math.jl")
include("utilities.jl")
include("state.jl")
include("nlsolver.jl")

# Models
include("models/structure.jl")
include("models/model.jl")
include("models/variational.jl")
include("models/single_polyion.jl")
include("models/coac_symmetric.jl")
include("models/coac_asymmetric.jl")
include("models/coac_association.jl")
include("models/coac_selfcomp.jl")
include("models/asymmetric.jl")

# Formulas
include("formulas/translational.jl")
include("formulas/chi.jl")
include("formulas/electrostatics.jl")
include("formulas/association.jl")
include("formulas/potentials.jl")
include("formulas/adaptive_chain.jl")

include("optimize.jl")
include("plotting.jl")

end