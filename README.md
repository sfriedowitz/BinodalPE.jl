# BinodalPE.jl

BinodalPE implements thermodynamic models (i.e. free energy functions) that allow for determination
of binodal coexistence conditions for polyelectrolyte solutions and coacervates.

Solving the binodal coexistence condition for a given model 
amounts to equating the electrochemical potentials of all species
and the osmotic pressure between the phases, overall mass balance by the lever rule, 
and charge neutrality in each macroscopic phase.

## Models

The package implements a series of `AbstractModel`s
that describe a specific physical system and free energy function.
Four main models are available:
1. `SinglePolyion` - A solution consisting of a single polyelectrolyte (A), counterions (+), and coions (-),
where the counterions can reversibly bind along the PE chain
2. `SymmetricCoacervate` - A solution consisting of a single "polymer" density (P) and a single "salt" density (S)
3. `AsymmetricCoacervate` - A solution with explicit treatment of polyanions (A), polycations (C), cations (+), and anions (-)
4. `AssociationCoacervate` - Same species as in an asymmetric coacervate, but with a thermodynamic description of local charge binding
This includes cation binding on the polyanions, anion binding along the polycations,
and interchain cross-linking betweeen the polyelectrolytes

Each model depends on numerous parameters that describe
the structural and energetic properties of salt and polyelectrolyte species,
including polymer chain lengths, charge densities, and monomer volumes.
These include: 
* `omega` - The normalized monomer volume of a species
* `sig` - The linear charge density of a polymer
* `chi` - The Flory-Huggins parameter between polymer and solvent
* `dp` - The degree of polymerization of a polymer
* `b` - The normalized segment length of a polymer
* `lp` - The reference persistence length of a semi-flexible polymer (only used for worm-like chains)
* `dg` - Short-range binding energy for models with reversible binding (`SinglePolyion` and `AssociationCoacervate`)

The electrostatic free energy within a model
depends on the conformational properties of the polyelectrolytes.
It is mandatory to specify an `AbstractChainStructure` when constructing a model
that specifies the structure factor of the polymers used within the electrostatic free energy.
Available chain structures include:
* `PointLike` - A Debye-Huckel description of point-like ions
* `ExtendedPoint` - The extended Debye-Huckel description for point-like ions
* `SmearedPoint` - A Debye-Huckel model for point-like ions with density smearing to regularize electrostatic interactions
* `GaussianCoil` - A flexible Gaussian coil described by the Debye function
* `RodLike` - A rigid, rod-like chain described by the Neugebauer function
* `WormLike` - A semiflexible, worm-like chain with an interpolated structure factor depending on the persistence length `lp`
* `SphericalGlobule` - A collapsed, globular polyelectrolyte
* `AdaptiveChain` - A worm-like chain with adaptively determined persistence length based on the RGF formalism (not available in all models)

## Usage

The general usage of this code is to:
1. Construct a model with a given set of parameters
2. Specify an initial bulk solution composition and initial guess for the compositions of coexisting phases
3. Sweep over a range of relevant parameters and solve for phase coexistence
with the `bndlsolve` or `bndlminimize` methods

This workflow is enabled by a series of general methods
implemented for each model, including 
`ftotal(phi, model)` which returns the free energy at a given composition `phi`,
`mutotal(phi, model)` which returns the chemical potential of each species evaluated at `phi`,
and `f2total(phi, model)` which returns the Hessian of the free energy evaluated at `phi`.

Additionally, for `SinglePolyion` and `AssociationCoacervate` models,
or any model with an `AdaptiveChain` structure specified,
one can solve for a set of internal state variables at a specified composition.
These may include binding fractions for reversible, short-range ion binding,
or the adaptively determined persistence length of an adaptive chain.
The internal solver routine for each model is implemented in the method `varsolve(phi, model)`.

Examples of these general procedures are shown below.

### Construct a model

```julia
julia> using BinodalPE 

julia> model_single = SinglePolyion(structure = GaussianCoil, omega = [5, 5, 1], dp = 100, sig = 1.0, dg = -5)

julia> model_symm = SymmetricCoacervate(structure = RodLike, omega = [5, 1], dp = 1000, sig = 0.5)

julia> modle_asymm = AsymmetricCoacervate(structure = WormLike, omega = [5, 5, 1, 1], dp = [100, 100], sig = [0.5, 0.5], lp = [10, 10])

julia> model_assoc = AssociationCoacervate(structure = AdaptiveChain, omega = [5, 5, 1, 1], dp = [100, 100], dg = [-4, -4, -8])
```

### Evaluate a model

```julia
julia> using BinodalPE

julia> model = SinglePolyion(structure = GaussianCoil, omega = [5, 1, 1], dp = 100, dg = -5)

julia> phi = [1e-4, 1.1e-3, 1e-3] # Volume fractions of polyanions, cations, and anions

julia> ftotal(phi, model) # Evaluate the solution free energy
-0.014242526701396035

julia> mutotal(phi, model) # Evaluate the chemical potentials of each species
3-element Array{Float64,1}:
 -1.0952518718803692
 -5.770961458699611
 -5.849320207862141

julia> varsolve(phi, model) # Solve for the binding fraction of counterions on the chain
1-element Array{Float64,1}:
 0.9241066860547466
```

### Solve a point on a phase diagram

```julia
julia> using BinodalPE

julia> model = SymmetricCoacervate(structure = GaussianCoil, dp = 100, sig = 0.25, omega = [1, 1])
SymmetricCoacervate{GaussianCoil}(bulk = [0.0, 0.0])

julia> phi = [0.01, 0.025]

julia> ftotal(phi, model)
-0.1017901482330647

julia> mutotal(phi, model)
2-element Array{Float64,1}:
 -0.6656002270678112
 -2.9963284857396184
 
julia> set_bulk!(model, [0.01, 0.072])

julia> init = [0.0005, 0.0725, 0.04056618971, 0.0725, 0.06] # Initial guesses for [phiP_sup, phiS_sup, phiP_coac, phiS_coac, volume_coac]
 
julia> res = bndlsolve(init, model; iterations = 50, ftol = 1e-10, show_trace = true)
Iter     f(x) inf-norm    Step 2-norm 
------   --------------   --------------
     0     2.137255e-02              NaN
     1     8.784988e-03     3.568122e+00
     2     2.164265e-03     1.042857e+00
     3     5.808519e-04     4.173605e-01
     4     5.729038e-05     1.050467e-01
     5     1.311902e-05     1.594545e-01
     6     5.933710e-07     3.427782e-02
     7     1.828610e-09     2.090991e-03
     8     1.489173e-14     5.844037e-06
BinodalResults:
  x = [0.006018056178356338, 0.07188967183864409, 0.017558300725491544, 0.07220941868076813, 0.34504847842390507]
  steps = 8
  objective = 1.4892e-14
  converged = true
  
julia> res.state
BinodalState:
  Bulk  = [0.01, 0.072]
  Sup   = [0.006018056178356338, 0.07188967183864409]
  Dense = [0.017558300725491544, 0.07220941868076813]
  ν = 0.34504847842390507
```

## TODO
* Expanded set of example notebooks for code usage and features
* Binodal Newton stability w/ sub optimal initial guesses
* Stable spinodal/critical point methods
* Transfer formulas into a type-based method for each interaction (modular interaction design)

## Publications

A number of publications describe the methods and theory implemented in this package.
We refer the user to these publications for details on theories implemented here.

1. Treatment of electrostatics
```bibtex
@article{Qin2016,
  doi = {10.1021/acs.macromol.6b02113},
  url = {https://doi.org/10.1021/acs.macromol.6b02113},
  year = {2016},
  month = nov,
  publisher = {American Chemical Society ({ACS})},
  volume = {49},
  number = {22},
  pages = {8789--8800},
  author = {Jian Qin and Juan J. de Pablo},
  title = {Criticality and Connectivity in Macromolecular Charge Complexation},
  journal = {Macromolecules}
}
```

2. Reversible charge association model
```bibtex
@article{Salehi2016,
  doi = {10.1021/acs.macromol.6b01464},
  url = {https://doi.org/10.1021/acs.macromol.6b01464},
  year = {2016},
  month = dec,
  publisher = {American Chemical Society ({ACS})},
  volume = {49},
  number = {24},
  pages = {9706--9719},
  author = {Ali Salehi and Ronald G. Larson},
  title = {A Molecular Thermodynamic Model of Complexation in Mixtures of Oppositely Charged Polyelectrolytes with Explicit Account of Charge Association/Dissociation},
  journal = {Macromolecules}
}
```
```bibtex
@article{Friedowitz2018,
  doi = {10.1063/1.5034454},
  url = {https://doi.org/10.1063/1.5034454},
  year = {2018},
  month = oct,
  publisher = {{AIP} Publishing},
  volume = {149},
  number = {16},
  pages = {163335},
  author = {Sean Friedowitz and Ali Salehi and Ronald G. Larson and Jian Qin},
  title = {Role of electrostatic correlations in polyelectrolyte charge association},
  journal = {The Journal of Chemical Physics}
}
```

3. Use of package to compute phase diagrams
```bibtex
@article{Lou2019,
  doi = {10.1021/acscentsci.8b00964},
  url = {https://doi.org/10.1021/acscentsci.8b00964},
  year = {2019},
  month = feb,
  publisher = {American Chemical Society ({ACS})},
  volume = {5},
  number = {3},
  pages = {549--557},
  author = {Junzhe Lou and Sean Friedowitz and Jian Qin and Yan Xia},
  title = {Tunable Coacervation of Well-Defined Homologous Polyanions and Polycations by Local Polarity},
  journal = {{ACS} Central Science}
}
```
