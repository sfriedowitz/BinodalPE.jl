@testset "state" begin

	state = BinodalState([0.01, 0.005])
	state.nu = 0.1

	new = copy(state)
	@test new.nu == state.nu

	cstate = toconc(state, [1, 1])
	@test cstate.bulk ≈ [0.5555555555555556, 0.2777777777777778]
	pstate = tophi(cstate, [1, 1])
	@test pstate.bulk ≈ [0.01, 0.005]

	data = BinodalData()
	for i = 1:10
		add_state!(data, BinodalState(rand(2)))
	end

	@test length(data) == 10
	deleteat!(data, 1)
	@test length(data) == 9
	
end