#==============================================================================#
# Single state
#==============================================================================#

"""
    mutable struct BinodalState

Holds information for a point on a binodal curve, including the `bulk` composition, 
a dilute `sup` state, a `dense` state, and the dense fraction `nu`.

Auxiliary properties can be added to the `props` dictionary field.
"""
mutable struct BinodalState{TB <: Real, TX <: Real, AB <: AbstractArray{TB}, AX <: AbstractArray{TX}}
    nu    :: TX
    bulk  :: AB
    sup   :: AX
    dense :: AX
    props :: Dict
end

BinodalState(bulk::AB) where AB = BinodalState(one(eltype(AB)), bulk, AB(zeros(length(bulk))), AB(zeros(length(bulk))), Dict())
BinodalState(bulk::AbstractArray{TB}, sup::AbstractArray{TX}, dense::AbstractArray{TX}, nu::TX) where {TB, TX} = BinodalState(nu, bulk, sup, dense, Dict())

function Base.show(io::IO, state::BinodalState)
    @printf(io, "BinodalState:\n")
    @printf(io, "  Bulk  = %s\n", state.bulk)
    @printf(io, "  Sup   = %s\n", state.sup)
    @printf(io, "  Dense = %s\n", state.dense)
    @printf(io, "  ν = %s", state.nu)
end

function Base.copy(state::BinodalState)
    new = BinodalState(copy(state.bulk), copy(state.sup), copy(state.dense), state.nu)
    for entry in state.props
        new.props[entry[1]] = entry[2]
    end
    return new
end

#==============================================================================#

function swap!(state::BinodalState)
    tmp = copy(state.sup)
    state.sup = state.dense
    state.dense = tmp
    state.nu = 1 - state.nu
    return nothing
end

function swap(state::BinodalState)
    new = copy(state)
    swap!(new)
    return new
end

function valid(state::BinodalState)
    # Check supernatant
    for v in state.sup
        if !inrange(v, 0.0, 1.0); return false; end
    end
    # Check dense phase
    for v in state.dense
        if !inrange(v, 0.0, 1.0); return false; end
    end
    # Check fraction
    if !inrange(state.nu, 0.0, 1.0)
        return false
    end
    # Check sums
    if sum(state.bulk) > 1.0 || sum(state.sup) > 1.0 || sum(state.dense) > 1.0
        return false
    end
    return true
end

function _convert_composition(state::BinodalState{TB,TX,AB,AX}, omega::AbstractArray{<:Real}, cmap) where {TB,TX,AB,AX}
    @assert length(omega) == length(state.bulk)

    new_bulk = convert(AB, cmap.(state.bulk, omega))
    new_sup = convert(AX, cmap.(state.sup, omega))
    new_dense = convert(AX, cmap.(state.dense, omega))

    new = BinodalState(new_bulk, new_sup, new_dense, state.nu)
    for (key, value) in state.props
        new.props[key] = value
    end

    return new
end

tophi(state::BinodalState, omega::AbstractArray{<:Real}) = _convert_composition(state, omega, tophi)
toconc(state::BinodalState, omega::AbstractArray{<:Real}) = _convert_composition(state, omega, toconc)

#==============================================================================#
# Holder for entire binodal
#==============================================================================#

"""
    struct BinodalData

Holder for multiple `BinodalState` points on a binodal curve.
"""
struct BinodalData
    states :: Vector{BinodalState}
    BinodalData() = new([])
end

Base.show(io::IO, bndl::BinodalData) = @printf(io, "BinodalData(%d states)", length(bndl.states))
Base.length(bndl::BinodalData) = length(bndl.states)
Base.iterate(bndl::BinodalData, idx = 1) = iterate(bndl.states, idx)
Base.reverse!(bndl::BinodalData) = reverse!(bndl.states)
Base.sort!(bndl::BinodalData, idx = 1) = sort!(bndl.states, by = x -> x.bulk[idx])
Base.empty!(bndl::BinodalData) = empty!(bndl.states)
Base.deleteat!(bndl::BinodalData, idx) = deleteat!(bndl.states, idx)
Base.getindex(bndl::BinodalData, idx::Integer) = getindex(bndl.states, idx)
Base.lastindex(bndl::BinodalData) = lastindex(bndl.states)

function Base.getindex(bndl::BinodalData, irange::OrdinalRange)
    new = BinodalData()
    for i in irange
        add_state!(new, bndl[i])
    end
    return new
end

#==============================================================================#

Base.push!(bndl::BinodalData, state::BinodalState) = push!(bndl.states, state)

function Base.push!(bndl::BinodalData, other::BinodalData)
    for state in other
        push!(bndl, state)
    end
    return nothing
end

add_state!(bndl::BinodalData, state::BinodalState) = push!(bndl.states, copy(state))

function tophi(bndl::BinodalData, omega::AbstractArray{<:Real})
    new = BinodalData()
    for s in bndl
        add_state!(new, tophi(s, omega))
    end
    return new
end

function toconc(bndl::BinodalData, omega::AbstractArray{<:Real})
    new = BinodalData()
    for s in bndl
        add_state!(new, toconc(s, omega))
    end
    return new
end

"""
    savebndl(bndl, fname)

Save a full binodal curve to JSON formatted file.
"""
function savebndl(bndl::BinodalData, fname)
    d = Dict()
    for (i, s) in enumerate(bndl)
        d[i] = Dict()
        
        d[i][:nu] = s.nu
        d[i][:bulk] = s.bulk
        d[i][:sup] = s.sup
        d[i][:dense] = s.dense
        d[i][:props] = s.props
    end

    json_string = JSON.json(d, 4)
    open(fname, "w") do f 
        write(f, json_string) 
    end

    return nothing
end

"""
    readbndl(bndl, fname)

Read a full binodal curve from a JSON formatted file.
"""
function readbndl(fname)
    d = JSON.parsefile(fname)
    dkeys = sort(parse.(Int, keys(d)))

    bndl = BinodalData()
    for k in dkeys
        vals = d[string(k)]
        
        nu = convert(Float64, get(vals, "nu", 0.0))
        bulk = convert(Vector{Float64}, get(vals, "bulk", []))
        sup = convert(Vector{Float64}, get(vals, "sup", []))
        dense = convert(Vector{Float64}, get(vals, "dense", []))
        props = vals["props"]
        
        state = BinodalState(bulk, sup, dense, nu)
        for entry in props
            if typeof(entry[2]) <: AbstractArray
                state.props[Symbol(entry[1])] = Float64.(entry[2])
            elseif typeof(entry[2]) <: Real
                state.props[Symbol(entry[1])] = Float64(entry[2])
            end
        end

        add_state!(bndl, state)
    end

    return bndl
end