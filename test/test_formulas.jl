@testset "formulas" begin

	# Test single polyion formulas
	phi = [0.01, 0.01, 0.01]
	model = SinglePolyion(structure = GaussianCoil, dg = -5)
	vars = varsolve(phi, model)

	@test vars[1] ≈ 0.8853291946844383
	@test ftotal(phi, model) ≈ -0.11991340426951363
	@test BinodalPE.felectrostatic(phi, vars, model) ≈ 0.012024772398863454

	# Test symmetric coacervate formulas
	phi = [0.01, 0.005]
	model = SymmetricCoacervate(structure = RodLike, omega = [5, 1])

	@test ftotal(phi, model) ≈ -0.03296564824229016
	@test BinodalPE.ftranslational(phi, model) ≈ -0.041470623479357406
	@test BinodalPE.felectrostatic(phi, model) ≈ 0.008504975237067242

	# Test asymmetric coacervate formulas
	phi = [0.01, 0.01, 0.005, 0.005]
	model = AsymmetricCoacervate(structure = WormLike, omega = [5, 5, 1, 1])

	@test ftotal(phi, model) ≈ -0.06530374287488896
	@test BinodalPE.ftranslational(phi, model) ≈ -0.0827128117330872
	@test BinodalPE.felectrostatic(phi, model) ≈ 0.017409068858198238

	# Test association coacervate formulas
	phi = [0.01, 0.01, 0.005, 0.005]
	model = AssociationCoacervate(structure = WormLike, omega = [5, 5, 1, 1])
	vars = varsolve(phi, model)

	@test vars ≈ [0.12293320886101733, 0.12293320886101725, 0.4330421273327021, 0.4330421273327021]
	@test ftotal(phi, model) ≈ -0.0679103840308165
	@test BinodalPE.ftranslational(phi, vars, model) ≈ -0.0805904282664085
	@test BinodalPE.felectrostatic(phi, vars, model) ≈ 0.012104535962343687
	
end