# SPDX-License-Identifier: MIT


function nikuradse(diameter, roughness)
    friction = (2 * log10(diameter / roughness) + 1.138)^(-2)
    @assert friction > 0
    return friction
end

function sigma(diameter, friction, c)
    σ = 16 * c^2 * friction / (π^2 * diameter^5)
    @assert σ > 0
    return σ
end

function theta_flux(diameter, friction)
    @assert diameter > 0 && friction > 0
    return friction / diameter
end

function get_model(solver, solvers_opts)
  logfile = nothing

  if solver == "Gurobi"
    model = direct_model(Gurobi.Optimizer())
    # Seems like a bug in Gurobi.jl, get_optimizer_attribute(m, "LogFile") fails
    # So use this kludge for now
    logfile = tempname()
    set_optimizer_attribute(model, "LogFile", logfile)
#  elseif solver == "CPLEX"
#    model = direct_model(CPLEX.Optimizer())
#  elseif solver == "Mosek"
#    # Mosek does not support ScalarQuadraticFunction in direct_model
#    model = Model(Mosek.Optimizer)
##    model = direct_model(Mosek.Optimizer())
  else
    msg = "ERROR: solver keyword \"$(solver)\" not supported"
    throw(ErrorException(msg))
  end

  # Transfer solver options
  if solver in keys(solvers_opts)
    for (k, v) in solvers_opts[solver]
      set_optimizer_attribute(model, k, v)
    end
  end

  return (model, logfile)
end
