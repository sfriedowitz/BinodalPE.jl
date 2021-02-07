@testset "binodal" begin

	# Test symmetric coacervate solver
	model = SymmetricCoacervate(structure = RodLike, sig = 0.25)
	set_bulk!(model, [0.01, 0.005])
	x0 = [1e-6, 0.004, 0.05, 0.006, 0.12]

	result = bndlsolve(x0, model, scale = true)
	@test result.x[end] ≈ 0.42625613931067097 atol = 1e-4

	# Test asymmetric coacervate solver
	model = AsymmetricCoacervate(structure = RodLike, sig = [0.25, 0.25])
	set_bulk!(model, [0.005, 0.005, 0.0025, 0.0025])
	x0 = [1e-8, 1e-8, 0.0025, 0.0025, 0.05, 0.05, 0.0025, 0.0025, 0.2, 0.0]

	result = bndlsolve(x0, model, scale = true)
	@test result.x[end-1] ≈ 0.42625613929281647 atol = 1e-4

	# Test association coacervate solver
	model = AssociationCoacervate(structure = GaussianCoil, dg = [-5, -5, -5], omega = [5, 5, 1, 1])
	set_bulk!(model, [0.005, 0.005, 0.0025, 0.0025])
	x0 = [1e-8, 1e-8, 0.0025, 0.0025, 0.05, 0.05, 0.0025, 0.0025, 0.2, 0.0]

	result = bndlsolve(x0, model, scale = true)
	@test result.x[end-1] ≈ 0.1267333610664382 atol = 1e-4
	
end