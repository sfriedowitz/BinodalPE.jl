# BinodalPE.jl

BinodalPE implements thermodynamic models (i.e. free energy functions) that allow for determination
of binodal coexistence conditions for polyelectrolyte solutions and coacervates.

Solving the binodal coexistence condition for a given model 
amounts to equating the electrochemical potentials of all species
and the osmotic pressure between the phases, overall mass balance by the lever rule, 
and charge neutrality in each macroscopic phase.

Four main models exist, each corresponding to a different physical system:
1. Single Polyion - A solution consisting of a single polyelectrolyte, counterions, and coions,
where the counterions can reversibly bind along the PE chain.
2. Symmetric Coacervate - A solution consisting of a single "polymer" density
and a single "salt" density, which can phase separate into a coexisting coacervate and supernatant.
3. Asymmetric Coacervate - A solution with explicit treatment of polyanions, polycations,
salt cations, and salt anions.
4. Association Coacervate - Same system as an asymmetric coacervate,
but where reversible charge association along the polyelectrolytes is treated.
Cation binding on the polyanions, anion binding along the polycations,
and interchain ion-pairing betweeen the PEs are considered.

An example of code usage is shown below:

```julia
julia> model = SymmetricCoacervate(structure = GaussianCoil, np = 100, sig = 0.25, omega = [1, 1])
SymmetricCoacervate{GaussianCoil}(bulk = [0.0, 0.0])

julia> phi = [0.01, 0.025];

julia> free_energy(phi, model)
-0.1017901482330647

julia> mutotal(phi, model)
2-element Array{Float64,1}:
 -0.6656002270678112
 -2.9963284857396184
 
julia> set_bulk(model, [0.01, 0.072])

julia> init = [0.0005, 0.0725, 0.04056618971, 0.0725, 0.06];
 
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
* Binodal Newton stability w/ sub optimal initial guesses
* Stable spinodal/critical point methods
* Transfer formulas into a type-based method for each interaction (modular interaction design)

## Publications

A number of publications describe the methods and theory implemented in this package.
We refer the user to these publications for details on the implementation.

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
