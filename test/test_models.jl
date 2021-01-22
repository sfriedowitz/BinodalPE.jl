@testset "models" begin

	models = [
		SinglePolyion(structure = GaussianCoil, dg = -5),
		SymmetricCoacervate(structure = RodLike, omega = [5, 1]),
		AsymmetricCoacervate(structure = WormLike, omega = [5, 5, 1, 1]),
		AssociationCoacervate(structure = SphericalGlobule, omega = [5, 5, 1, 1], dg = [-5, -5, -5])
	]

	phis = [
	    [0.01, 0.01, 0.01],
	    [0.01, 0.005],
	    [0.01, 0.01, 0.005, 0.005],
	    [0.01, 0.01, 0.005, 0.005]
	]

	for (i, model) in enumerate(models)
		@test newstate(model).nu == 1.0

		set_bulk(model, phis[i])
		@test model.bulk == phis[i]
	end
	
end