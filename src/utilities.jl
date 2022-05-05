#==============================================================================#
# Helper functions and macros
#==============================================================================#

notimpl(name, type) = error("Error: `$(name)` not implemented for $(type).")

inrange(x::T, a = T(0.0), b = T(1.0)) where T <: Number = (a < x <= b)

function array_perturb!(x::AbstractArray, delta::Real = 0.0)
	for i in eachindex(x)
		x[i] += delta * randn()
	end
	return nothing
end

function logscale!(x::AbstractArray, irange::AbstractUnitRange = eachindex(x))
	for i in irange
		x[i] = log(x[i]/(1 - x[i]))
	end
	return nothing
end

logscale(x::AbstractArray, irange::AbstractUnitRange = eachindex(x)) = (new = copy(x); logscale!(new, irange); new)

logscale(x::StaticArray) = @. log(x/(1-x))

function logunscale!(x::AbstractArray, irange::AbstractUnitRange = eachindex(x))
	for i in irange
		x[i] = 1/(1 + exp(-x[i]))
	end
	return nothing
end

logunscale(x::AbstractArray, irange::AbstractUnitRange = eachindex(x)) = (new = copy(x); logunscale!(new, irange); new)

logunscale(x::StaticArray) = @. 1/(1 + exp(-x))

#==============================================================================#

toconc(p, w, v = v0) = p / (1000.0*Navo*w*v)

tophi(c, w, v = v0) = 1000.0 * (Navo*c*w*v)

asypolyion(pb, r) = (pb*(1+r)/2, pb*(1-r)/2)