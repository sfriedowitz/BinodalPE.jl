using Test
using BinodalPE

#==============================================================================#

add_jl(x) = endswith(x, ".jl") ? x : x*".jl"

if length(ARGS) > 0
    tests = map(add_jl, ARGS)
else
    tests = [
        "test_models.jl",
        "test_formulas.jl",
        "test_binodal.jl",
        "test_state.jl"
    ]
end

println("Testing...")

for test in tests
    include(test)
end