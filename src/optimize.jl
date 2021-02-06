#==============================================================================#
# Results storage
#==============================================================================#

struct BinodalResults
    x         :: Vector{Float64}
    bulk      :: Vector{Float64}
    state     :: BinodalState
    steps     :: Int
    objective :: Float64
    converged :: Bool
end

function Base.show(io::IO, res::BinodalResults)
    @printf(io, "BinodalResults:\n")
    @printf(io, "  x = %s\n", res.x)
    @printf(io, "  steps = %d\n", res.steps)
    @printf(io, "  objective = %.5g\n", res.objective)
    @printf(io, "  converged = %s\n", res.converged)
end

#==============================================================================#
# Generic Newton scheme for solving coexistence condition
#==============================================================================#

function bndlf(x, model::AbstractModel; scaled::Bool = false)
    F = similar(x)
    bndlf!(F, x, model; scaled = scaled)
    return F
end

function bndlj(x, model::AbstractModel; scaled::Bool = false)
    J = zeros(eltype(x), length(x), length(x))
    bndlj!(J, x, model; scaled = scaled)
    return J
end

function bndlfj!(F, J, x, model::AbstractModel; scaled::Bool = false)
    bndlf!(F, x, model; scaled = scaled)
    bndlj!(J, x, model; scaled = scaled)
    return nothing
end

#==============================================================================#

bndlsolve(res::BinodalResults, model::AbstractModel; kwargs...) = bndlsolve(res.x, model; kwargs...)
bndlsolve(state::BinodalState, model::AbstractModel; kwargs...) = bndlsolve(bndlsolvex(state, model), model; kwargs...)

"""
    bndlsolve(init, model; kwargs...)
    bndlsolve(results, model; kwargs...)
    bndlsolve(state, model; kwargs...)

Solve for binodal coexistence by equating chemical potentials and pressures between each phase.
Initial parameters provided in `init` are in an order specific to each `model`.

Keword arguments:
* `tf`             : Floating-point type used in calculations
* `iterations`     : Maximum number of Newton iterations (50)
* `xtol`           : Stopping tolerance for change in x vector (0.0)
* `ftol`           : Stopping tolerance for residual vector norm (1e-10)
* `rlxn`           : Relaxation factor for Newton solver (0.9)
* `pmax`           : Maximum Newton step magnitude (Inf)
* `scale`          : Perform optimization in scaled variable format (true)
* `show_trace`     : Print iteration history (false)
* `extended_trace` : Print extended trace if `show_trace` is true (false)
* `linesearch`     : Linesearch method from `LineSearches.jl` to use in solver (BackTracking(order = 3))
"""
function bndlsolve(init::AbstractVector, model::AbstractModel;
    tf::Type{TF} = Float64, 
    iterations::Integer = 50,
    xtol::Real = 0.0,
    ftol::Real = 1e-10,
    rlxn::Real = 0.9,
    pmax::Real = Inf,
    scale::Bool = true,
    show_trace::Bool = false,
    extended_trace::Bool = false,
    autodiff::Symbol = :auto,
    linesearch = LineSearches.Static()
) where {TF <: Real}

    # Initial vectors
    x0 = convert(Vector{TF}, copy(init))
    if scale; bndlscale!(x0, model); end
    F0 = similar(x0)

    f!(F, x) = bndlf!(F, x, model; scaled = scale)
    j!(J, x) = bndlj!(J, x, model; scaled = scale)
    fj!(F, J, x) = bndlfj!(F, J, x, model; scaled = scale)
    if autodiff == :finite
        df = OnceDifferentiable(f!, x0, F0, autodiff = autodiff)
    else
        df = OnceDifferentiable(f!, j!, fj!, x0, F0)
    end

    # Internal solver routine after setup
    sol = _newton_solve(df, x0, iterations, xtol, ftol, rlxn, pmax, false, show_trace, extended_trace, linesearch)
    #sol = nlsolve(df, x0, method = :trust_region, iterations = iterations, ftol = ftol, xtol = xtol, show_trace = show_trace, extended_trace = extended_trace, linesearch = linesearch)

    # State and results
    xsol = bndlunscale(copy(sol.zero), model)
    return BinodalResults(xsol, copy(model.bulk), bndlstate(xsol, model), sol.iterations, sol.residual_norm, sol.f_converged)
end

#==============================================================================#
# Free energy minimization method
#==============================================================================#

bndlminimize(res::BinodalResults, model::AbstractModel; kwargs...) = bndlminimize(res.x, model; kwargs...)
bndlminimize(state::BinodalState, model::AbstractModel; kwargs...) = bndlminimize(bndlminx(state, model), model; kwargs...)

"""
    bndlminimize(init, model; kwargs...)
    bndlminimize(results, model; kwargs...)
    bndlminimize(state, model; kwargs...)

Perform a flash calculation to obtain the compositions of coexisting phases 
at the bulk concentration specified in the `model`. The order of variables in the initial guess `init` 
is unique to each model.

An optimization routine from the `Optim.jl` package is provided in the argument `method`.

Keword arguments:
* `method`     : Optimization method to perform (NelderMead())
* `cycles`     : Number of minimization cycles (1)
* `iterations` : Number of steps per cycle (500)
* `delta`      : Parameter perturbation scaling between cycles (1e-5)
* `perturb`    : Turn on perturbation between cycles (true)
* `show_trace` : Print trace of minimization cycles (false)
"""
function bndlminimize(init::AbstractVector, model::AbstractModel;
    method::Optim.AbstractOptimizer = NelderMead(),
    cycles::Integer = 1,
    iterations::Integer = 500,
    ftol::Real = 1e-10,
    delta::Real = 1e-5,
    perturb::Bool = true,
    show_trace::Bool = false
)

    # Initial vectors
    x = copy(init)
    state = bndlstate(x, model)
    if !valid(state)
        error("Error: Initial variables produce invalid state!")
    end
    bndlscale!(x, model)

    # Setup objective function
    fb = free_energy(model.bulk, model)
    fx = x -> bndlg(bndlunscale(x, model), model; fbulk = fb)
    options = Optim.Options(f_tol = ftol, iterations = iterations)

    xmin = copy(x)
    fmin = fx(xmin)
    cycle_min = 0

    if show_trace
        @printf("Starting minimization -- f(x) = %.5g\n", fmin)
    end

    itx = 0
    converged = false
    for c = 1:cycles
        res = optimize(fx, x, method, options)
        itx += res.iterations

        if res.minimum < fmin
            fmin = res.minimum
            x = copy(res.minimizer)
            xmin = copy(x)
            converged = res.f_converged

            cycle_min = c
        end

        if show_trace
            @printf("Minimization cycle (%d)  -- f(x) = %.5g\n", c, fmin)
        end

        state = bndlstate(bndlunscale(x, model), model)
        if state.sup[1] > state.dense[1]
            if show_trace
                println("Swapping dilute and dense phases...")
            end

            swap!(state)
            x = bndlminx(state, model)
            bndlscale!(x, model)

            if cycle_min == c
                xmin = copy(x)
            end
        end

        if perturb && c < cycles
            array_perturb!(x, delta)
        end
    end

    sol = bndlunscale(xmin, model)
    state = bndlstate(sol, model)
    return BinodalResults(sol, copy(model.bulk), state, itx, fmin, converged)
end