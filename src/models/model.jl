#==============================================================================#
# Abstract model type
#==============================================================================#

"""
    abstract type AbstractModel{TC <: AbstractChainStructure}

A thermodynamic model of a molecular system.
"""
abstract type AbstractModel{TC <: AbstractChainStructure} end

Base.show(io::IO, model::TM) where TM <: AbstractModel = @printf(io, "%s(bulk = %s)", TM, model.bulk)

#==============================================================================#

"""
    neutralbulk(phi_bulk, model)

Return an appropriate bulk composition with the added salt and polymer values in `phi_bulk`,
accounting for charge neutrality and the presence of counterions from the polymers.
"""
neutralbulk(phi, model::AbstractModel) = notimpl("neutralbulk", typeof(model))

function set_bulk!(model::AbstractModel, bulk)
    # Some error checking
    if length(bulk) != length(model.bulk)
        error("Error: Invalid number of bulk components, requires $(length(model.bulk)).")
    elseif sum(bulk) > 1.0
        error("Error: Bulk composition exceeds unity.")
    end
    model.bulk .= bulk
    return nothing
end

#==============================================================================#

"""
    bndlscale(x, model)

Return scaled coexistence variables for an `AbstractModel`
to improve numerical stability for binodal solvers.
"""
function bndlscale(x, model::AbstractModel)
    xs = copy(x)
    bndlscale!(xs, model)
    return xs
end

bndlscale!(x, model::AbstractModel) = notimpl("bndlscale!", typeof(model))

"""
    bndlunscale(x, model)

Return unscaled coexistence variables for an `AbstractModel` to original form.
"""
function bndlunscale(xs, model::AbstractModel)
    x = copy(xs)
    bndlunscale!(x, model)
    return x
end

bndlunscale!(x, model::AbstractModel) = notimpl("bndlunscale!", typeof(model))

#==============================================================================#

newstate(model::AbstractModel) = BinodalState(model.bulk)

chains(model::AbstractModel) = notimpl("chains(model)", typeof(model))

chains(model::AbstractModel, vars) = chains(model)

bndlminx(state::BinodalState, model::AbstractModel) = notimpl("bndlminx", typeof(model))

bndlsolvex(state::BinodalState, model::AbstractModel) = notimpl("bndlsolvex", typeof(model))

bndlstate(x, model::AbstractModel) = notimpl("bndlstate", typeof(model))

bndlf!(F, x, model::AbstractModel) = notimpl("bndlf!", typeof(model))

bndlj!(J, x, model::AbstractModel) = notimpl("bndlj!", typeof(model))

function bndlg(x, model::AbstractModel; fbulk::Real = NaN)
    state = bndlstate(x, model)
    fb = !isnan(fbulk) ? fbulk : ftotal(state.bulk, model)

    if !valid(state)
        return 1e10*rand()
    else
        fs = ftotal(state.sup, model)
        fc = ftotal(state.dense, model)
        return ((1-state.nu)*fs + state.nu*fc - fb)/abs(fb)
    end
end

bndlg(state::BinodalState, model::AbstractModel) = bndlg(bndlminx(state, model), model)

#==============================================================================#
# Free energy and derivatives
#==============================================================================#

"""
    ftotal(phi, model)

Return the total solution free energy for a model at a given composition.
"""
ftotal(phi, model::AbstractModel) = notimpl("ftotal", typeof(model))

ftranslational(phi, model::AbstractModel) = notimpl("ftranslational", typeof(model))

function fideal(phi, model::AbstractModel)
    phiW = 1-sum(phi)
    return ftranslational(phi, model) - phiW*log(phiW)
end

function fexcess(phi, model::AbstractModel)
    phiW = 1-sum(phi)
    ftot = ftotal(phi, model)
    return ftot - ftranslational(phi, model) + phiW*log(phiW)
end

function f2total(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    # Without scaling procedure
    gx = x -> ftotal(x, model)
    cfg = ForwardDiff.HessianConfig(gx, phi, chunk)
    return ForwardDiff.hessian(gx, phi, cfg)
end

function f3total(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    n = length(phi)
    fx = x -> f2total(x, model; chunk = chunk)
    cfg = ForwardDiff.JacobianConfig(fx, phi, chunk)
    J = ForwardDiff.jacobian(fx, phi, cfg)
    return reshape(J, (n,n,n))
end

function mutotal(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    # # Without scaling procedure
    # gx = x -> ftotal(x, model)
    # cfg = ForwardDiff.GradientConfig(gx, phi, chunk)
    # return ForwardDiff.gradient(gx, phi, cfg)

    # With variable scaling
    xs = logscale(phi)
    fx = x -> ftotal(logunscale(x), model)

    cfg = ForwardDiff.GradientConfig(fx, xs, chunk)
    g = ForwardDiff.gradient(fx, xs, cfg)

    dx = @. 1/(phi - phi^2)
    return g .* dx
end

function muideal(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    fx = x -> fideal(x, model)
    cfg = ForwardDiff.GradientConfig(fx, phi, chunk)
    return ForwardDiff.gradient(fx, phi, cfg)
end

function muexcess(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    fx = x -> fexcess(x, model)
    cfg = ForwardDiff.GradientConfig(fx, phi, chunk)
    return ForwardDiff.gradient(fx, phi, cfg)
end

function pressure(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    mu = mutotal(phi, model; chunk = chunk)
    return dot(mu, phi) - ftotal(phi, model)
end

function mupressure(phi, model::AbstractModel; chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(phi))
    mu = mutotal(phi, model; chunk = chunk)
    p = dot(mu, phi) - ftotal(phi, model)
    return (mu, p)
end