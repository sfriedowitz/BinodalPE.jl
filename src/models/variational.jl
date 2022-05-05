#==============================================================================#
# Variational solver functions
#==============================================================================#

"""
    varinit(phi, model)

Provide an initial guess for an `AbstractModel`'s internal self-consistent solver routine.
"""
varinit(phi, model::AbstractModel) = notimpl("varinit", typeof(model))

"""
    varscale(x, model)

Rescale variational parameters, by default according to x -> log(x/(1-x))
"""
varscale(x, model::AbstractModel) = logscale(x)

"""
    varunscale(x, model)

Rescale variational parameters, by default according to x -> exp(x)/(1 + exp(x))
"""
varunscale(x, model::AbstractModel) = logunscale(x)

#==============================================================================#

"""
    varf!(F, x, phi, model)

Fill the objective function `F` for the variational method defined for a `model`.
"""
varf!(F, x, phi, model::AbstractModel) = notimpl("varf!", typeof(model))

"""
    varf(x, phi, model; scale = false)

Return the objective function `F` for the `model`'s variational solver
given the state variables `x`.
"""
function varf(x, phi, model::AbstractModel; scale::Bool = false)
    F = similar(x)
    xs = scale ? varscale(x, model) : x
    varf!(F, xs, phi, model)
    return F
end

function varj!(J, x, phi, model::AbstractModel)
    F = similar(x)
    fx! = (F, z) -> varf!(F, z, phi, model)
    ForwardDiff.jacobian!(J, fx!, F, x)
    return nothing
end

function varfj!(F, J, x, phi, model::AbstractModel)
    fx! = (F, z) -> varf!(F, z, phi, model)
    ForwardDiff.jacobian!(J, fx!, F, x)
    return nothing
end

#==============================================================================#

"""
    varsolve(phi, model)

Self-consistently solve for the internal state variables for an `AbstractModel` 
at a given composition `phi`.
"""
function varsolve(phi, model::AbstractModel;
    iterations::Integer = 50,
    xtol::Real = 0.0,
    ftol::Real = 1e-10,
    rlxn::Real = 1.0,
    pmax::Real = 100.0,
    autodiff::Symbol = :finite,
    scaling::AbstractArray = [],
    store_trace::Bool = false,
    show_trace::Bool = false,
    extended_trace::Bool = false,
    linesearch = LineSearches.BackTracking(order = 3),
    linsolve = (x, A, b) -> copyto!(x, A\b),
)
    # Setup vectors and df
    init = varinit(phi, model)
    x0 = varscale(init, model)
    F0 = similar(x0)
    TF = eltype(x0)

    function f!(F, x)
        varf!(F, x, phi, model)
        for i = 1:min(length(F), length(scaling))
            F[i] *= scaling[i]
        end
        return nothing
    end
    df = OnceDifferentiable(f!, x0, F0, autodiff)

    # Internal solver routine after setup
    sol = _newton_solve(df, x0, iterations, convert(TF, xtol), convert(TF, ftol), convert(TF, rlxn), convert(TF, pmax),
                        store_trace, show_trace, extended_trace, linesearch, linsolve)

    # check for convergence
    (sol.residual_norm <= ftol && !isnan(sol.residual_norm)) || println("varsolve not converged for phi = $phi with $(typeof(model)), f(x) = $(sol.residual_norm) ")

    return varunscale(sol.zero, model)
end