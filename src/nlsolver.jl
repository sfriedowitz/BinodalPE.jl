#==============================================================================#
# Modified Newton algorithm
#==============================================================================#

import NLsolve: OnceDifferentiable, NewtonCache, SolverTrace, SolverResults
import NLsolve: newtontrace, check_isfinite, assess_convergence
import NLsolve: value, jacobian, value!, value_jacobian!, value_jacobian!!

function _newton_solve(df::OnceDifferentiable, initial_x::AbstractArray{T},
    iterations::Integer,
    xtol::Real,
    ftol::Real,
    rlxn::Real,
    pmax::Real,
    store_trace::Bool,
    show_trace::Bool,
    extended_trace::Bool,
    linesearch = LineSearches.BackTracking(order = 2),
    linsolve = (x, A, b) -> copyto!(x, A\b),
    cache = NewtonCache(df)
) where T

    # Setup cache and vectors
    n = length(initial_x)
    copyto!(cache.x, initial_x)
    value_jacobian!!(df, cache.x)
    check_isfinite(value(df))
    vecvalue = vec(value(df))
    x_ls = copy(cache.x)

    it = 0
    stopped = any(isnan, cache.x) || any(isnan, value(df)) ? true : false
    x_converged, f_converged = assess_convergence(cache.x, cache.xold, value(df), xtol, ftol)
    converged = x_converged || f_converged

    # Setup tracing
    if show_trace
        @printf "Iter     f(x) inf-norm    Step 2-norm \n"
        @printf "------   --------------   --------------\n"
    end
    tr = SolverTrace()
    tracing = store_trace || show_trace || extended_trace
    newtontrace(convert(real(T), NaN), tracing, extended_trace, cache, df, it, tr, store_trace, show_trace)

    # Create objective function for the linesearch.
    # This function is defined as fo(x) = 0.5 * f(x) ⋅ f(x) and thus
    # has the gradient ∇fo(x) = ∇f(x)' ⋅ f(x)
    # where ∇f(x) is the Jacobian matrix of the problem
    function fo!(a)
        @. cache.x = cache.xold + a*cache.p
        value!(df, cache.x)
        return dot(value(df), value(df))/2
    end
    function go!(a)
        @. cache.x = cache.xold + a*cache.p
        value_jacobian!(df, cache.x)
        mul!(cache.g, transpose(jacobian(df)), vecvalue)
        return dot(cache.g, cache.p)
    end
    function fgo!(a)
        @. cache.x = cache.xold + a*cache.p
        value_jacobian!(df, cache.x)
        mul!(cache.g, transpose(jacobian(df)), vecvalue)
        return dot(value(df), value(df))/2, dot(cache.g, cache.p)
    end
    
    while !stopped && !converged && it < iterations
        it += 1
        if it > 1
            # Linear mixing of previous states -- seems to improve stability
            @. cache.x = (1 - rlxn)*cache.xold + rlxn*cache.x

            # Update F and J
            value_jacobian!(df, cache.x)
        end

        try
            mul!(vec(cache.g), transpose(jacobian(df)), vec(value(df)))
            linsolve(cache.p, jacobian(df), vec(value(df)))
            rmul!(cache.p, -1)
        catch e
            if isa(e, LAPACKException) || isa(e, SingularException)
                # Modify the search direction if the jacobian is singular
                # FIXME: better selection for lambda, see Nocedal & Wright p. 289
                J2 = transpose(jacobian(df))*jacobian(df)
                lam = convert(real(T),1e6)*sqrt(n*eps())*norm(J2, 1)
                linsolve(cache.p, -(J2 + lam * I), value(df))
            else
                throw(e)
            end
        end

        # Set limit of step size magnitude before linesearch
        if norm(cache.p) > pmax
            cache.p ./= norm(cache.p)
            cache.p .*= pmax
        end

        # Backup old state
        copyto!(cache.xold, cache.x)

        # Perform linesearch
        if linesearch isa Static
            alpha = one(real(T))
            phi_alpha = fo!(alpha)
        else
            alpha, phi_alpha  = linesearch(fo!, go!, fgo!, one(real(T)), dot(value(df), value(df))/2, dot(cache.g, cache.p))
        end

        # x and F updated in linesearch, so no need to call again
        stopped = any(isnan, cache.x) || any(isnan, value(df)) ? true : false
        x_converged, f_converged = assess_convergence(cache.x, cache.xold, value(df), xtol, ftol)
        converged = x_converged || f_converged

        newtontrace(sqeuclidean(cache.x, cache.xold), tracing, extended_trace, cache, df, it, tr, store_trace, show_trace)
    end

    return SolverResults("Newton with line-search",
                         initial_x, copy(cache.x), norm(value(df), Inf),
                         it, x_converged, xtol, f_converged, ftol, tr,
                         first(df.f_calls), first(df.df_calls))
end