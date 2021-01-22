# """
#     bndlminimize(init, model; kwargs...)
#     bndlminimize(results, model; kwargs...)
#     bndlminimize(state, model; kwargs...)

# Perform a flash calculation to obtain the compositions of coexisting phases 
# at the bulk concentration specified in the `model`.

# Keword arguments:
# * `optimizer`  : An instance of `MOI.OptimizerWithAttributes` to solve the problem. `Ipopt.Optimizer` is recommended for nonlinear problems.
# * `cycles`     : Repeat optimization cycles (1)
# * `delta`      : Parameter perturbation scaling between cycles (1e-5)
# * `extra_cons` : Add extra bounding constraints to the problem (true)
# """
# function bndlminimize(init::AbstractVector, model::AbstractModel;
#     optimizer::MOI.OptimizerWithAttributes,
#     cycles::Integer = 1,
#     delta::Real = 1e-5,
#     extra_cons::Bool = true,
#     show_trace::Bool = false
# )
#     # Create JuMP nonlinear modeling problem
#     problem = bndlproblem(model; init = bndlscale(init, model), extra_constraints = extra_cons)
#     set_optimizer(problem, optimizer)

#     # Run cycles of optimization
#     fsol = Inf
#     psol = copy(init)

#     # Setup tracing
#     if show_trace
#         @printf "Nonlinear Optimization Cycles\n"
#     end

#     for c in 1:cycles
#         JuMP.optimize!(problem)

#         # Take if better
#         if objective_value(problem) < fsol
#             psol = JuMP.value.(problem[:p])
#             fsol = objective_value(problem)
#         end

#         if show_trace
#             @printf "  Cycle %d: f(x) = %.5g\n" c fsol
#         end

#         # Perturb starting points
#         if cycles > 1 && c < cycles
#             reinit = bndlscale(psol, model)
#             perturb!(reinit, delta)
#             set_start_value.(problem[:x], reinit)
#         end
#     end

#     converged = termination_status(problem) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
#     return BinodalResults(psol, copy(model.bulk), bndlstate(psol, model), cycles, fsol, converged)
# end

#==============================================================================#

# function bndlproblem(model::SymmetricCoacervate; init::Union{AbstractArray,Nothing} = nothing, extra_constraints::Bool = true)
#     # Generate model
#     problem = JuMP.Model()

#     # Parameters
#     @NLparameter(problem, b[i=1:2] == model.bulk[i])

#     # Variables
#     @variable(problem, x[i=1:5], start = x0[i])
#     @NLexpression(problem, p[i=1:5], exp(x[i])/(1 + exp(x[i])))

#     if !isnothing(init)
#         for i in 1:length(x)
#             set_start_value(x[i], init[i])
#         end
#     end

#     # Add constraints
#     @NLconstraint(problem, [i=1:2], b[i] == (1-p[5])*p[i] + p[5]*p[i+2])
#     if extra_constraints
#         @NLconstraint(problem, [i=1:5], 0 <= p[i] <= 1)
#         @NLconstraint(problem, p[3] >= p[1])
#     end
    
#     # Objective function
#     fb = free_energy(model.bulk, model)
#     gx(p...) = bndlg(p, model; fbulk = fb)
#     JuMP.register(problem, :gx, 5, gx, autodiff = true)

#     @NLobjective(problem, Min, gx(p...))
    
#     return problem
# end

# function bndlproblem(model::AsymmetricCoacervateModel; init::Union{AbstractArray,Nothing} = nothing, extra_constraints::Bool = true)
#     # Generate model
#     problem = JuMP.Model()

#     # Parameters    
#     @NLparameter(problem, b[i=1:4] == model.bulk[i])
#     @NLparameter(problem, w[i=1:4] == model.omega[i])
#     @NLparameter(problem, f[i=1:2] == model.sig[i])

#     # Variables
#     @variable(problem, x[i=1:9])
#     @NLexpression(problem, p[i=1:9], exp(x[i])/(1 + exp(x[i])))

#     if !isnothing(init)
#         for i in 1:length(x)
#             set_start_value(x[i], init[i])
#         end
#     end

#     # Add constraints
#     @NLconstraint(problem, [i=1:4], b[i] == (1-p[9])*p[i] + p[9]*p[i+4])
#     @NLconstraint(problem, -f[1]*p[1]/w[1] + f[2]*p[2]/w[2] + p[3]/w[3] - p[4]/w[4] == 0)
#     @NLconstraint(problem, -f[1]*p[5]/w[1] + f[2]*p[6]/w[2] + p[7]/w[3] - p[8]/w[4] == 0)
    
#     if extra_constraints
#         @NLconstraint(problem, [i=1:9], 0 <= p[i] <= 1)
#         @NLconstraint(problem, p[1] <= p[5])
#         @NLconstraint(problem, p[2] <= p[6])
#     end

#     # Objective function
#     fb = free_energy(model.bulk, model)
#     gx(p...) = bndlg(p, model; fbulk = fb)
#     JuMP.register(problem, :gx, 9, gx, autodiff = true)

#     @NLobjective(problem, Min, gx(p...))
    
#     return problem
# end

# function update_state!(x, state::BinodalState, model::AsymmetricCoacervateModel)
#     if length(x) == 4
#         # Offload parameters
#         phiAC, phiCC, phiWC, nu = x
#         phiAB, phiCB, phiPB, phiMB = model.bulk
#         phiWB = 1 - sum(model.bulk)
#         wA, wC, wP, wM = model.omega
#         sigA, sigC = model.sig
        
#         # Derive dense phase params
#         phiPC = wP/(wP+wM) - phiAC*(wP/(wP+wM)) - phiCC*(wP/(wP+wM)) - phiWC*(wP/(wP+wM)) + sigA*phiAC*(wM*wP/(wA*(wP+wM))) - sigC*phiCC*(wM*wP/(wC*(wP+wM)))
#         phiMC = 1 - phiAC - phiCC - phiPC - phiWC
        
#         # Derive supernatant phase params
#         phiAS = (phiAB - nu*phiAC)/(1 - nu)
#         phiCS = (phiCB - nu*phiCC)/(1 - nu)
#         phiWS = (phiWB - nu*phiWC)/(1 - nu)
#         phiPS = wP/(wP+wM) - phiAS*(wP/(wP+wM)) - phiCS*(wP/(wP+wM)) - phiWS*(wP/(wP+wM)) + sigA*phiAS*(wM*wP/(wA*(wP+wM))) - sigC*phiCS*(wM*wP/(wC*(wP+wM)))
#         phiMS = 1 - phiAS - phiCS - phiPS - phiWS

#         # Update state struct
#         state.nu = nu
#         state.sup .= (phiAS, phiCS, phiPS, phiMS)
#         state.dense .= (phiAC, phiCC, phiPC, phiMC)

#     elseif length(x) == 10
#         # Offload parameters
#         phiAB, phiCB, phiPB, phiMB = model.bulk
#         phiAS, phiCS, phiPS, phiMS = x[1], x[2], x[3], x[4]
#         phiAC, phiCC, phiPC, phiMC = x[5], x[6], x[7], x[8]
#         nu, psi = x[9], x[10]

#         # Update state struct
#         state.nu = nu
#         state.bulk .= (phiAB, phiCB, phiPB, phiMB)
#         state.sup .= (phiAS, phiCS, phiPS, phiMS)
#         state.dense .= (phiAC, phiCC, phiPC, phiMC)
#         state.props[:psi] = psi
#     end
#     return nothing
# end