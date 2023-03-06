# SPDX-License-Identifier: MIT
using Crayons.Box, Dates, HDF5, JSON, JuMP, LinearAlgebra, Logging

include("helpers.jl")
include("readers.jl")
include("utils.jl")
include("writers.jl")
include("logger.jl")
include("econ_computations.jl")

g = 9.81

macro central_diff_2nd_order(var, idx)
  return :($var[$idx+1] + $var[$idx-1] - 2*$var[$idx])
end

macro fwd_diff_1st_order(var, idx)
  return :($var[$idx+1] - $var[$idx])
end


function solve(m, params)

  area = params["area"]
  c = params["speed_of_sound"]
  θ = params["theta_flux"]
  g = params["gravity"]

  pipe_data = params["pipe_data"]
  slope = pipe_data["slope"]

  nb_t_pts = params["nb_t_pts"]
  nb_x_pts = params["nb_x_pts"]
  Δ = params["delta_s"]
  Δt = params["delta_t"]
  Δx = params["delta_x"]

  econ_data = params["econ_data"]
  FricL = params["FricL"]

  pressure_lb = params["pressure_lb"]
  pressure_ub = params["pressure_ub"]
  flux_like_lb = params["flux_like_lb"]
  flux_like_ub = params["flux_like_ub"]

  scale_R = params["scale_R"]

  model_opts = params["model_opts"]
  penalization_choices = get(model_opts, "penalization_choices", [])

  task_durations = Dict{String,Float64}()

  if isempty(penalization_choices)
    penalization_choice = get(model_opts, "penalization_choice", nothing)
    if penalization_choice !== nothing
      penalization_choices = [penalization_choice]
    end
  end


  terminal_condition_choice = get(model_opts, "terminal_condition_choice", "none")
  control_boundary = get(model_opts, "control_boundary", "hyperbolic_compatible")

  smoothing_pts = get(model_opts, "smoothing_pts", 0)
  smoothing_delta = Δt * smoothing_pts

  row_scaling = Dict{String,Float64}()

  # In Julia, range(a, length=len) returns [0, len-1]
  time_idxs_bndry = range(0, length=nb_t_pts)
  time_idxs_full = range(0, length=2*nb_t_pts-1)
  x_idxs = range(0, length=nb_x_pts)

#  x_pts = range(0, params["pipe_data"]["length"], nb_x_pts)
#

  ###############################
  # Step 1: define the variables
  ###############################

  @info GREEN_FG("Defining variables")
  t = time()

  # fully_free is the parabolic case
  if control_boundary == "fully_free"
    ctrl_idxs = (x_idxs[1], x_idxs[end])
    @variable(m, Rp_ctrl[time_idxs_bndry, ctrl_idxs] >= 0)
    @variable(m, Rm_ctrl[time_idxs_bndry, ctrl_idxs] <= 0)
  else
    @variable(m, Rp_ctrl[time_idxs_bndry] >= 0)
    @variable(m, Rm_ctrl[time_idxs_bndry] <= 0)
  end

  @variable(m, Rp[time_idxs_full, x_idxs] >= 0)
  @variable(m, Rm[time_idxs_full, x_idxs] <= 0)

  task_duration = round(time() - t, digits=1)
  task_durations["def:vars"] = task_duration
  @info GREEN_FG("Done in $(task_duration) secs")

  #####################################
  # Step 2: define all the constraints
  #####################################


  @info GREEN_FG("Defining the constraints")
  t = time()
  # Impose initial conditions
  Rp_initial, Rm_initial = params["Riemann_init"]

  # Note that the "+1" is here since Julia is 1-based
  @constraint(m, Rp_init[i=x_idxs], Rp[0, i] == Rp_initial[i+1])
  @constraint(m, Rm_init[i=x_idxs], Rm[0, i] == Rm_initial[i+1])

  # Impose the controls at each end of the pipe

  if control_boundary == "fully_free"
    @constraint(m, Rp_control[t_idx=time_idxs_bndry[2:end],x_idx=ctrl_idxs], Rp[2*t_idx, x_idx] == Rp_ctrl[t_idx,x_idx])
    @constraint(m, Rm_control[t_idx=time_idxs_bndry[2:end],x_idx=ctrl_idxs], Rm[2*t_idx, x_idx] == Rm_ctrl[t_idx,x_idx])

    # At t = 0, we impose that the control should be compatible with the initial conditions
    @constraint(m, Rp_control_init[x_idx=ctrl_idxs], Rp_ctrl[0,x_idx] == Rp_initial[1+x_idx])
    @constraint(m, Rm_control_init[x_idx=ctrl_idxs], Rm_ctrl[0,x_idx] == Rm_initial[1+x_idx])
  else
    @constraint(m, Rp_control[t_idx=time_idxs_bndry[2:end]], Rp[2*t_idx, 0] == Rp_ctrl[t_idx])
    @constraint(m, Rm_control[t_idx=time_idxs_bndry[2:end]], Rm[2*t_idx, end] == Rm_ctrl[t_idx])

    # At t = 0, we impose that the control should be compatible with the initial conditions
    @constraint(m, Rp_control_init, Rp_ctrl[0] == Rp_initial[1])
    @constraint(m, Rm_control_init, Rm_ctrl[0] == Rm_initial[end])
  end

  # Define the different relations depending on the even or odd value of the time index

  # Rp is fixed at x = 0 
  rge_Rp_x_even  = range(1, length=nb_x_pts-1)
  # Rm is fixed at x = L 
  rge_Rm_x_even  = range(0, length=nb_x_pts-1)

  # For odd lines, we actually have one less variable
  rge_Rp_x_odd   = range(0, length=nb_x_pts-1)
  rge_Rm_x_odd   = range(0, length=nb_x_pts-1)

  # For the time dynamic, it's way easier
  rge_t_idx_odd  = range(1, step=2, length=nb_t_pts-1)
  rge_t_idx_even = range(2, step=2, length=nb_t_pts-1)

  # Fix unnecessary variables to a feasible value
  fix.(Rp[rge_t_idx_odd, end], (pressure_lb + pressure_ub)/2 + (flux_like_lb + flux_like_ub)/4; force=true)
  fix.(Rm[rge_t_idx_odd, end], (flux_like_lb + flux_like_ub)/4 - (pressure_lb + pressure_ub)/2; force=true)

  #############################################
  # Add the relations representing the dynamics
  #############################################

  # We scale the relation to fit into the range recommended by Gurobi, see
  # https://www.gurobi.com/documentation/9.5/refman/advanced_user_scaling.html
  # Our target is to get the coefficient in the range [1e-3, 1e6]
  #
  coeff_dyn = 1/sqrt(Δ/2 *  θ*c/4)
  row_scaling["_dyn_"] = coeff_dyn


  # Define Rp for even time index
  # Rpk = Rp[t_idx, x_idx]
  # Rmk = Rm[t_idx, x_idx]
  # We have the convention that odd lines have their last component fixed
  # Rpi = Rp[t_idx-1, x_idx-1]
  # Rmi = Rm[t_idx-1, x_idx-1]
  @constraint(m, Rp_dyn_even[t_idx=rge_t_idx_even,x_idx=rge_Rp_x_even],
                 coeff_dyn * (Rp[t_idx, x_idx] - Rp[t_idx-1, x_idx-1]) == 
                 coeff_dyn * (-Δ/2)*( θ*c/4  * (FricL(m, Rp, Rm, t_idx, x_idx) + FricL(m, Rp, Rm, t_idx-1, x_idx-1))
                      + g * slope / (2*c) * (Rp[t_idx, x_idx] - Rm[t_idx, x_idx] + Rp[t_idx-1, x_idx-1] - Rm[t_idx-1, x_idx-1])))

  # Define Rm for even time index
  # Rpk = Rp[t_idx, x_idx]
  # Rmk = Rm[t_idx, x_idx]
  # We have the convention that odd lines have their last component fixed
  # Rpj = Rp[t_idx-1, x_idx]
  # Rmj = Rm[t_idx-1, x_idx]
  @constraint(m, Rm_dyn_even[t_idx=rge_t_idx_even,x_idx=rge_Rm_x_even],
                 coeff_dyn * (Rm[t_idx, x_idx] - Rm[t_idx-1, x_idx]) == 
                 coeff_dyn * (-Δ/2)*( θ*c/4  * (FricL(m, Rp, Rm, t_idx, x_idx) + FricL(m, Rp, Rm, t_idx-1, x_idx))
                      + g * slope / (2*c) * (Rp[t_idx, x_idx] - Rm[t_idx, x_idx] + Rp[t_idx-1, x_idx] - Rm[t_idx-1, x_idx])))


  # Define Rp for odd time index
  # Rpk = Rp[t_idx, x_idx]
  # Rmk = Rm[t_idx, x_idx]
  # Rpi = Rp[t_idx-1, x_idx]
  # Rmi = Rm[t_idx-1, x_idx]
  @constraint(m, Rp_dyn_odd[t_idx=rge_t_idx_odd,x_idx=rge_Rp_x_odd],
                 coeff_dyn * (Rp[t_idx, x_idx] - Rp[t_idx-1, x_idx]) ==
                 coeff_dyn * (-Δ/2)*( θ*c/4  * (FricL(m, Rp, Rm, t_idx, x_idx) + FricL(m, Rp, Rm, t_idx-1, x_idx))
                      + g * slope / (2*c) * (Rp[t_idx, x_idx] - Rm[t_idx, x_idx] + Rp[t_idx-1, x_idx] - Rm[t_idx-1, x_idx])))

  # Define Rm for odd time index
  # Rpk = Rp[t_idx, x_idx]
  # Rmk = Rm[t_idx, x_idx]
  # Rpj = Rp[t_idx-1, x_idx+1]
  # Rmj = Rm[t_idx-1, x_idx+1]
  @constraint(m, Rm_dyn_odd[t_idx=rge_t_idx_odd,x_idx=rge_Rm_x_odd],
                 coeff_dyn * (Rm[t_idx, x_idx] - Rm[t_idx-1, x_idx+1]) == 
                 coeff_dyn * (-Δ/2)*( θ*c/4  * (FricL(m, Rp, Rm, t_idx, x_idx) + FricL(m, Rp, Rm, t_idx-1, x_idx+1))
                         + g * slope / (2*c) * (Rp[t_idx, x_idx] - Rm[t_idx, x_idx] + Rp[t_idx-1, x_idx+1] - Rm[t_idx-1, x_idx+1])))


  # CPLEX does not support interval constraints
  rge_Rp_tidxs, rge_Rp_xidxs = axes(Rp)
  @constraint(m, pressure_lb[t_idx=rge_Rp_tidxs,x_idx=rge_Rp_xidxs], pressure_lb <= (Rp[t_idx, x_idx] - Rm[t_idx, x_idx])/2)
  @constraint(m, pressure_ub[t_idx=rge_Rp_tidxs,x_idx=rge_Rp_xidxs], (Rp[t_idx, x_idx] - Rm[t_idx, x_idx])/2 <= pressure_ub)
  @constraint(m, mflow_lb[t_idx=rge_Rp_tidxs,x_idx=rge_Rp_xidxs], flux_like_lb <= Rp[t_idx, x_idx] + Rm[t_idx, x_idx])
  @constraint(m, mflow_ub[t_idx=rge_Rp_tidxs,x_idx=rge_Rp_xidxs], Rp[t_idx, x_idx] + Rm[t_idx, x_idx] <= flux_like_ub)

  rge_t_idx_even_middle = range(2, step=2, length=nb_t_pts-2)
  # common factor: transformation from flux to mass flow
  # To try to provide a better scaling, we do not multiply by the constant
  # cst = Δt * 2 * c * area
  @constraint(m, preservation_gas_rule,
              sum(Rp[tk, 0] + Rm[tk, 0] for tk in rge_t_idx_even_middle) + (Rp[0, 0] + Rm[0, 0] + Rp[end, 0] + Rm[end, 0])/2 >=
              sum(Rp[tk, end] + Rm[tk, end] for tk in rge_t_idx_even_middle) + (Rp[0, end] + Rm[0, end] + Rp[end, end] + Rm[end, end])/2)
  row_scaling["preservation_gas_rule"] =  Δt * 2 * c * area



  booking_number = get(model_opts, "booking_constraint", NaN)
  if isfinite(booking_number)
    cst = Δt * area / (2*c * scale_R)
    @constraint(m, booking, sum(Rp_ctrl[:, end] + Rm_ctrl[:, end])
                            - .5*( Rp_ctrl[0, end] + Rm_ctrl[0, end]
                                 + Rp_ctrl[end, end] + Rm_ctrl[end, end]) <= booking_number/cst)
  end

  terminal_penalization_control_expr = @expression(m, 0)

  if (terminal_condition_choice == "terminal_mflow_identical")

    @constraint(m, terminal_mflow_identical[x_idx=range(0, length=nb_x_pts-1)],
                Rp[end, x_idx] + Rm[end, x_idx] == Rp[end, end] + Rm[end, end])

  elseif (terminal_condition_choice == "terminal_curvature_penalization")

    penalty_terminal_curvature = params["num_coeffs"]["penalty_terminal_curvature"]

    terminal_penalization_control_expr = @expression(m, penalty_terminal_curvature / Δx^4 * (sum(
        ((Rp[end,idx+1] + Rm[end,idx+1]) + (Rp[end,idx-1] + Rm[end,idx-1]) - 2*(Rp[end,idx] + Rm[end,idx]))^2
        for idx in range(1, stop=nb_x_pts-2))))

  elseif (terminal_condition_choice == "terminal_gradient_penalization")

    penalty_terminal_gradient = params["num_coeffs"]["penalty_terminal_gradient"]

    terminal_penalization_control_expr = @expression(m, penalty_terminal_gradient / Δx^2 * (sum(
        ((Rp[end,idx+1] + Rm[end,idx+1]) - (Rp[end,idx] + Rm[end,idx]))^2
        for idx in range(0, stop=nb_x_pts-2))))

  elseif (terminal_condition_choice == "none")
    # all fine here
  else
    println("ERROR: unsupported terminal_condition_choice \"$(terminal_condition_choice)\"")
    exit()
  end

  task_duration = round(time() - t, digits=1)
  task_durations["def:cons"] = task_duration
  @info GREEN_FG("Done in $(task_duration) secs")

  ########################################
  # Step 3: define the objective function
  ########################################

  @info GREEN_FG("Defining objective function")
  t = time()

  Δt = params["delta_t"]
  time_pts = range(econ_data["time_intervals"][1], params["T_c"], length=nb_t_pts)
  opt_mflows = Vector{Float64}(undef, nb_t_pts-1)
  obj_econ_exprs = Vector{Any}(undef, nb_t_pts-1)

  for (idx, tk) in enumerate(time_pts[1:end-1])
    tkp1 = time_pts[idx + 1]

    @assert abs(tkp1 - tk - Δt) < 1e-10

    vidx = idx-1

#    revenue_lin_term, revenue_sqr_term, production_costs_term = get_obj_csts(econ_data, tk, tkp1)
    econ_csts = get_obj_csts_v2(econ_data, tk, tkp1; smoothing_delta=smoothing_delta)
    revenue_cst_tk   = econ_csts["revenue_cst_tk"]
    revenue_cst_tkp1 = econ_csts["revenue_cst_tkp1"]
    revenue_lin_tk   = econ_csts["revenue_lin_tk"]
    revenue_lin_tkp1 = econ_csts["revenue_lin_tkp1"]
    prodcost_cst_tk   = econ_csts["prodcost_cst_tk"]
    prodcost_cst_tkp1 = econ_csts["prodcost_cst_tkp1"]
    prodcost_lin_tk   = econ_csts["prodcost_lin_tk"]
    prodcost_lin_tkp1 = econ_csts["prodcost_lin_tkp1"]

    # Convert reduced flux (scaled by scale_R and in kg/(s * m^2)) to mass flow (kg/s)
    mflow_in_tk    = area*(Rp[2*vidx, 0] + Rm[2*vidx, 0])/(2*c * scale_R)
    mflow_in_tkp1  = area*(Rp[2*(vidx+1), 0] + Rm[2*(vidx+1), 0])/(2*c * scale_R)
    mflow_out_tk   = area*(Rp[2*vidx, end] + Rm[2*vidx, end])/(2*c * scale_R)
    mflow_out_tkp1 = area*(Rp[2*(vidx+1), end] + Rm[2*(vidx+1), end])/(2*c * scale_R)

    # Note that the data from get_obj_csts_v2 is already multiplied by period_len
    obj_econ_exprs[idx] = @expression(m, 1/2 * ( (revenue_cst_tk + revenue_lin_tk * mflow_out_tk ) * mflow_out_tk
                                               + (revenue_cst_tkp1 + revenue_lin_tkp1 * mflow_out_tkp1 ) * mflow_out_tkp1
                                               - (prodcost_cst_tk + prodcost_lin_tk * mflow_in_tk) * mflow_in_tk
                                               - (prodcost_cst_tkp1 + prodcost_lin_tkp1 * mflow_in_tkp1) * mflow_in_tkp1) )


    # Here we approximate somehow further by taking the averages of the coefficient
    revenue_cst_avg = (revenue_cst_tk + revenue_cst_tkp1)/2 
    revenue_lin_avg = (revenue_lin_tk + revenue_lin_tkp1)/2 
    prodcost_cst_avg = (prodcost_cst_tk + prodcost_cst_tkp1)/2 
    prodcost_lin_avg = (prodcost_lin_tk + prodcost_lin_tkp1)/2 
    opt_mflows[idx] = (prodcost_cst_avg - revenue_cst_avg)/(2*(revenue_lin_avg-prodcost_lin_avg))
    #println("Time interval [$(tk), $(tkp1)]: optimal stationary mass flow $(opt_mflows[idx])")

  end

  penalization_control_expr_laplacian = @expression(m, 0)
  penalization_control_expr_grad = @expression(m, 0)
  penalization_control_expr_L2 = @expression(m, 0)
  # It seems that the use of macro in a sum is a no go :(
  if "symmetric difference" in penalization_choices
    penality_scaling = params["num_coeffs"]["penalty_obj_ctrl_laplacian"]
    idxs = time_idxs_bndry[2:end-1]
    partial_exprs = Vector{Any}(undef, length(idxs))
    lidx = 1
    if control_boundary == "fully_free"
      for idx in idxs
        partial_exprs[lidx] = @expression(m, (Rp_ctrl[idx+1, ctrl_idxs[1]] + Rp_ctrl[idx-1, ctrl_idxs[1]] - 2*Rp_ctrl[idx, ctrl_idxs[1]])^2
                                           + (Rp_ctrl[idx+1, ctrl_idxs[2]] + Rp_ctrl[idx-1, ctrl_idxs[2]] - 2*Rp_ctrl[idx, ctrl_idxs[2]])^2
                                           + (Rm_ctrl[idx+1, ctrl_idxs[1]] + Rm_ctrl[idx-1, ctrl_idxs[1]] - 2*Rm_ctrl[idx, ctrl_idxs[1]])^2
                                           + (Rm_ctrl[idx+1, ctrl_idxs[2]] + Rm_ctrl[idx-1, ctrl_idxs[2]] - 2*Rm_ctrl[idx, ctrl_idxs[2]])^2)
      end
    else
      for idx in idxs
        partial_exprs[lidx] = @expression(m, (Rp_ctrl[idx+1] + Rp_ctrl[idx-1] - 2*Rp_ctrl[idx])^2 + (Rm_ctrl[idx+1] + Rm_ctrl[idx-1] - 2*Rm_ctrl[idx])^2)
        lidx += 1
      end
    end
    penalization_control_expr_laplacian = @expression(m, penality_scaling * Δt/Δt^4 * sum(ex for ex in partial_exprs))
    println("Magnitude for symmetric difference is $(penality_scaling * Δt/Δt^4)")
  end

  if "forward difference" in penalization_choices
    penality_scaling = params["num_coeffs"]["penalty_obj_ctrl_grad"]
    println("Magnitude for forward difference is $(penality_scaling * Δt/Δt^2)")
    idxs = time_idxs_bndry[1:end-1]
    if control_boundary == "fully_free"
      penalization_control_expr_grad = @expression(m, penality_scaling * Δt/Δt^2 *
        ( sum((Rp_ctrl[idx+1, ctrl_idxs[1]] - Rp_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rp_ctrl[idx+1, ctrl_idxs[2]] - Rp_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs)
        + sum((Rm_ctrl[idx+1, ctrl_idxs[1]] - Rm_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rm_ctrl[idx+1, ctrl_idxs[2]] - Rm_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs) ))

    else
      penalization_control_expr_grad = @expression(m, penality_scaling * Δt/Δt^2 * (sum((Rp_ctrl[idx+1] - Rp_ctrl[idx])^2 for idx in idxs)
                                                                      + sum((Rm_ctrl[idx+1] - Rm_ctrl[idx])^2 for idx in idxs)))
    end
  end

  if "forward difference (p,q)" in penalization_choices
    penality_scaling_p = params["num_coeffs"]["penalty_obj_ctrl_grad_p"]
    penality_scaling_q = params["num_coeffs"]["penalty_obj_ctrl_grad_q"]
    coeff = Δt/(Δt^2 * scale_R)
    println("Magnitude for forward difference for p is $(penality_scaling_p/2 * coeff) and q is $(area/(2*c) * penality_scaling_q* coeff)")
    idxs = time_idxs_bndry[1:end-1]
    if control_boundary == "fully_free"
      penalization_control_expr_grad =
      @expression(m, Δt/(Δt^2 * scale_R) * (area/(2*c) * penality_scaling_q *
        ( sum((Rp_ctrl[idx+1, ctrl_idxs[1]] + Rm_ctrl[idx+1, ctrl_idxs[1]] -
               Rp_ctrl[idx, ctrl_idxs[1]] - Rm_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rp_ctrl[idx+1, ctrl_idxs[2]] + Rm_ctrl[idx+1, ctrl_idxs[2]] -
               Rp_ctrl[idx, ctrl_idxs[2]] - Rm_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs)))
      + ( penality_scaling_p / 2 *
        ( sum((Rp_ctrl[idx+1, ctrl_idxs[1]] - Rm_ctrl[idx+1, ctrl_idxs[1]] -
               Rp_ctrl[idx, ctrl_idxs[1]] + Rm_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rp_ctrl[idx+1, ctrl_idxs[2]] - Rm_ctrl[idx+1, ctrl_idxs[2]] -
               Rp_ctrl[idx, ctrl_idxs[2]] + Rm_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs))
       ))

    else
      penalization_control_expr_grad = @expression(m, penality_scaling * Δt/Δt^2 * (sum((Rp_ctrl[idx+1] - Rp_ctrl[idx])^2 for idx in idxs)
                                                                      + sum((Rm_ctrl[idx+1] - Rm_ctrl[idx])^2 for idx in idxs)))
    end
  end

  if "L2" in penalization_choices
    # This would be an implicit approximation to get things simple
    penality_scaling = params["num_coeffs"]["L2"]
    println("Magnitude for forward difference is $(penality_scaling * Δt)")
    idxs = time_idxs_bndry[1:end]
    if control_boundary == "fully_free"
      penalization_control_expr_L2 = @expression(m, penality_scaling * Δt *
        ( sum((Rp_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rp_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs)
        + sum((Rm_ctrl[idx, ctrl_idxs[1]])^2 for idx in idxs)
        + sum((Rm_ctrl[idx, ctrl_idxs[2]])^2 for idx in idxs) ))

    else
      penalization_control_expr_L2 = @expression(m, penality_scaling * Δt *
        ( sum((Rp_ctrl[idx])^2 for idx in idxs)
        + sum((Rm_ctrl[idx])^2 for idx in idxs) ))
    end
  end

  valid_pen = ["symmetric difference", "forward difference", "forward difference (p,q)", "L2"]
  for pen in penalization_choices
    if pen ∉ valid_pen
      msg = "Unknown penalization_choice value $pen. "
      msg += "It should be one of the followng $(join(valid_pen, ", " ))"
      throw(ErrorException(msg))
    end
  end

  obj_scaling_factor = params["num_coeffs"]["scaling_obj"]

  penalization_control_expr = penalization_control_expr_L2 + penalization_control_expr_grad + penalization_control_expr_laplacian

  profit_expr = @expression(m, sum(ex for ex in obj_econ_exprs))
  @objective(m, Max, obj_scaling_factor * (profit_expr - penalization_control_expr - terminal_penalization_control_expr))

  task_duration = round(time() - t, digits=1)
  task_durations["def:obj"] = task_duration
  @info GREEN_FG("Done in $(task_duration) secs")


  @info GREEN_FG("Optimizing")
  t = time()

  optimize!(m)

  task_duration = round(time() - t, digits=1)
  task_durations["optimize"] = task_duration
  @info GREEN_FG("Done in $(task_duration) secs")

  return Dict{String,Any}("opt_mflows" => opt_mflows,
                          "row_scaling" => row_scaling,
                          "profit_expr" => profit_expr,
                          "penalization_control_expr" => penalization_control_expr,
                          "terminal_penalization_control_expr" => terminal_penalization_control_expr)
end


function save_data(m::JuMP.Model, fid::HDF5.File, params::Dict{String,Any};
                   save_multipliers::Bool=false, row_scaling::Dict{String,Float64}=Dict{String,Float64}())

  num_coeffs = params["num_coeffs"]
  c = params["speed_of_sound"]
  area = params["area"]
  Δt = params["delta_t"]
  scale_R = params["scale_R"]

  ##############################################################
  # Step 1: First get variables and their domain from the model
  ##############################################################

  Rp = m[:Rp]
  Rm = m[:Rm]
  Rp_ctrl = m[:Rp_ctrl]
  Rm_ctrl = m[:Rm_ctrl]
  time_idxs_bndry, = axes(Rp_ctrl)
  time_idxs_full, x_idxs = axes(Rp)

  RpIN = value.(Rp[2*time_idxs_bndry, 0]).data
  RmIN = value.(Rm[2*time_idxs_bndry, 0]).data

  RpOUT = value.(Rp[2*time_idxs_bndry, end]).data
  RmOUT = value.(Rm[2*time_idxs_bndry, end]).data

  h5bndry = create_group(fid, "/boundary")
  h5ctrl = create_group(fid, "/control")
  h5pipe = create_group(fid, "/pipe")
  h5units = create_group(fid, "/units")

  h5units["pressure"] = "bar"
  h5units["Riemann variable"] = "bar"
  h5units["mflow"] = "kg / s"

  h5bndry["RmIN"] = RmIN
  h5bndry["RpIN"] = RpIN
  h5bndry["RmOUT"] = RmOUT
  h5bndry["RpOUT"] = RpOUT

  pIN, fluxIN = Riemann2physics(RpIN, RmIN, c)
  pOUT, fluxOUT = Riemann2physics(RpOUT, RmOUT, c)

  h5ctrl["pIN"] = pIN
  h5ctrl["mflowIN"] = area*fluxIN/scale_R
  h5ctrl["pOUT"] = pOUT
  h5ctrl["mflowOUT"] = area*fluxOUT/scale_R

  p_min = Inf
  p_max = -Inf
  q_min = Inf
  q_max = -Inf

  Rp_mult_min = Inf
  Rp_mult_max = -Inf
  Rm_mult_min = Inf
  Rm_mult_max = -Inf

  if save_multipliers
    h5ctrl_mult = create_group(fid, "/control_mult")
    h5bnds_mult = create_group(fid, "/bounds_mult")

    h5ctrl_mult["Rp"] = dual.(m[:Rp_control]).data
    h5ctrl_mult["Rm"] = dual.(m[:Rm_control]).data

    # Scale the Lagrange multipliers according to the row scaling
    scaled_cons = keys(row_scaling)
    scale_pressure_bnd = 1.
    scale_dyn = 1.
    for (k,v) in row_scaling
      if occursin(k, "pressure_")
        scale_pressure_bnd = v
      elseif occursin(k, "_dyn_even")
        scale_dyn = v
      end
    end

    h5bnds_mult["pLB"] = scale_pressure_bnd * dual.(m[:pressure_lb]).data
    h5bnds_mult["pUB"] = scale_pressure_bnd * dual.(m[:pressure_ub]).data
    h5bnds_mult["mflowLB"] = dual.(m[:mflow_lb]).data
    h5bnds_mult["mflowUB"] = dual.(m[:mflow_ub]).data


    Rp_dyn_even = m[:Rp_dyn_even]
    Rm_dyn_even = m[:Rm_dyn_even]
  end

  for idx in time_idxs_bndry

    h5pipe_tk = create_group(h5pipe, string(idx))
    h5pipe_tk["time"] = idx * Δt

    full_idx = 2*idx
    Rp_pipe = value.(Rp[full_idx, :]).data
    Rm_pipe = value.(Rm[full_idx, :]).data
    # TODO add option to disable this to save space
    h5pipe_tk["Rp"] = Rp_pipe
    h5pipe_tk["Rm"] = Rm_pipe

    p_pipe_bar, flux_pipe_scaled = Riemann2physics(Rp_pipe, Rm_pipe, c)

    h5pipe_tk["pressure"] = p_pipe_bar
    h5pipe_tk["mflow"] = mflow = area*flux_pipe_scaled/scale_R

    # recompute the extremas
    p_min_tk, p_max_tk = extrema(p_pipe_bar)
    p_min = min(p_min, p_min_tk)
    p_max = max(p_max, p_max_tk)

    q_min_tk, q_max_tk = extrema(mflow)
    q_min = min(q_min, q_min_tk)
    q_max = max(q_max, q_max_tk)

    if save_multipliers && full_idx >= 2
      Rp_mult = scale_dyn * dual.(Rp_dyn_even[full_idx, :]).data
      Rm_mult = scale_dyn * dual.(Rm_dyn_even[full_idx, :]).data
      h5pipe_tk["Rp_mult"] = Rp_mult
      h5pipe_tk["Rm_mult"] = Rm_mult

      Rp_mult_min_tk, Rp_mult_max_tk = extrema(Rp_mult)
      Rp_mult_min = min(Rp_mult_min, Rp_mult_min_tk)
      Rp_mult_max = max(Rp_mult_max, Rp_mult_max_tk)
      Rm_mult_min_tk, Rm_mult_max_tk = extrema(Rm_mult)
      Rm_mult_min = min(Rm_mult_min, Rm_mult_min_tk)
      Rm_mult_max = max(Rm_mult_max, Rm_mult_max_tk)
    end

  end

  # Save top-level values

  fid["pressure_range"] = [p_min, p_max]
  fid["mflow_range"] = [q_min, q_max]

  fid["Rp_mult_range"] = [Rp_mult_min, Rp_mult_max]
  fid["Rm_mult_range"] = [Rm_mult_min, Rm_mult_max]


  fid["solvestat"] = Int(termination_status(m))
  fid["solve_time"] = solve_time(m)
  fid["barrier_iter"] = barrier_iterations(m)
  fid["simplex_iter"] = simplex_iterations(m)

  fid["obj_value"] = objective_value(m) / num_coeffs["scaling_obj"]

  close(fid)
end



function hPDE_DO(;cfg_dir::String="../cfg/example-debug",
                  nb_x_pts::Int64=50,
                  tmp_dir::String="/tmp",
                  json_data::Dict{String,Any}=Dict{String,Any}(),
                  solvers_opts::Dict{String,Any}=Dict{String,Any}(),
                  model_opts::Dict{String,Any}=Dict{String,Any}(),
                  solver="Gurobi",
                  output_dir::String="",
                  save_multipliers::Bool=true)

  model, logfile = get_model(solver, solvers_opts)

  # Allow to take the data from arguments
  if all(map(x -> haskey(json_data, x), ["pipe.json", "econ.json", "approx.json"]))
    pipe_data = json_data["pipe.json"]
    econ_data = json_data["econ.json"]
    approx_data = json_data["approx.json"]
  else
    pipe_data = pipe_reader(joinpath(cfg_dir, "pipe.json"))
    econ_data = econ_reader(joinpath(cfg_dir, "econ.json"))
    approx_data = approx_reader(joinpath(cfg_dir, "approx.json"))
  end

  # Open the HDF5 file
  timestamp = Dates.format(Dates.now(), "yyyy-mm-ddTHH_MM_SS")
  mkpath(output_dir)
  hdf5_filename = output_dir * "result-$(timestamp).hdf5"
  fid = hdf5_writer_init(hdf5_filename, pipe_data, econ_data, approx_data)

  logio = IOBuffer()
  global_logger(get_logger(logio))

  ##################################
  # End of the initialization part #
  ##################################
  
  params = Dict{String, Any}("model_opts" => model_opts)

  params["pipe_data"] = pipe_data
  params["econ_data"] = econ_data
  params["approx_data"] = approx_data

  L = params["pipe_data"]["length"]

  speed_of_sound = sqrt(pipe_data["Rs"] * pipe_data["T"])
  params["speed_of_sound"] = c = speed_of_sound
  params["gravity"] = g

  params["nb_x_pts"] = nb_x_pts
  params["delta_x"] = Δx = L / (nb_x_pts - 1)
  # this is the time difference between an even and odd row
  params["delta_s"] = Δx / (2*c)
  D = pipe_data["diameter"]
  params["area"] = area = π * (D/2)^2

  # time difference between consecutive nodes with same spatial position
  params["delta_t"] = Δt = Δx / c

  # Number of nodes on the I_mid segment
  @assert params["econ_data"]["time_intervals"][1] == 0.

  T = params["econ_data"]["time_intervals"][end]
  I_mid = T - L/c

  @assert I_mid > 0.

  # N_mid includes the nodes at the 2 vertices
  N_mid = ceil(Int64, I_mid/Δt) + 1

  # We readjust the total time to make sure that the boundaries are all aligned
  params["T_c"] = T_c = L/c + (N_mid-1) * Δt

  if T_c < T
    error("T_c = $(T_c) < $(T) = T")
  end

  # odd but correct
  params["nb_t_pts"] = nb_t_pts = N_mid + nb_x_pts - 1


  num_coeffs_default = get_num_coeffs_default(Δx, Δt)
  num_coeffs = get(model_opts, "num_coeffs", nothing)
  if isnothing(num_coeffs)
    num_coeffs = num_coeffs_default
  else
    for (k, v) in num_coeffs_default
      if k ∉ keys(num_coeffs)
        num_coeffs[k] = v
      end
    end
  end

  diff = setdiff(keys(num_coeffs), keys(num_coeffs_default))
  if length(diff) > 0
    @error "unknown numerical coefficient named $(join(diff, ", "))"
    exit()
  end


  # Scaling
  params["scale_R"] = scale_R = 1e-5

  # pOP and qOP are SI
  qOP = approx_data["operational_massflow"]
  pOP = approx_data["operational_pressure"]
  # flux = mflow/area
  fluxOP = qOP/area
 
  friction = nikuradse(D, pipe_data["k"])
  params["theta_flux"] = theta_flux(D, friction)

  params["pressure_lb"] = scale_R * pipe_data["pressure_min"]
  params["pressure_ub"] = scale_R * pipe_data["pressure_max"]

  # flux = mflow/area
  # Then q_flux * 2c = (R_+)  +  (R_-)
  params["flux_like_lb"] = scale_R * 2*c * pipe_data["flow_min"] / area
  params["flux_like_ub"] = scale_R * 2*c * pipe_data["flow_max"] / area


  σ = sigma(D, friction, c)
  # This is Weymouth
  # p(x) ^2 - pin^2 = - Λ q|q| x
  init_pressure_Pa(idx) = (pOP^2 - qOP * abs(qOP) * σ * idx * Δx)^(.5)

  init_stationary_state = get(model_opts, "init_stationary_equation", "LISO2")
  if init_stationary_state == "LISO2"
    pressure_init = scale_R*p_init_LISO2(pOP, qOP, c, collect(range(0,length=nb_x_pts) * Δx), pipe_data)
  elseif init_stationary_state == "ISO2"
    pressure_init = scale_R*init_pressure_Pa.(range(0,length=nb_x_pts))
  else
    @error "incorrect \"init_stationary_equation\" parameter value in model_opts: $(init_stationary_state) ∉ [LISO2, ISO2]"
    exit()
  end


  flux_init = scale_R*fluxOP*ones(nb_x_pts)
  params["FricL"] = get_linearized_friction_term(pOP, fluxOP, c, scale_R)

  Rp_init, Rm_init = physic2Riemann(pressure_init, flux_init, c)

  params["Riemann_init"] = (Rp_init, Rm_init)
  params["num_coeffs"] = num_coeffs

  seq_iter = 1
  seq_iter_max = get(model_opts, "seq_iter_max", 1)
  seq_err_tol = get(model_opts, "seq_err_tol", 1e-2)
  RpHat_old = nothing
  RmHat_old = nothing
  RpHat_oldold = nothing
  RmHat_oldold = nothing
  p_oldold = nothing
  q_oldold = nothing
  solve_artifacts = nothing
  norm_err = get(num_coeffs, "normErr", Inf)

  # Start logging
  grp_log = create_group(fid, "/log")
  grp_loop = hdf5_init_loop_log(grp_log, "PDE linearization"; tol=seq_err_tol, norm=norm_err, iter_max=seq_iter_max)

  errp = NaN
  errm = NaN
  err_pressure = NaN
  err_flux = NaN

  while true
    solve_artifacts = solve(model, params)
    RpHat = value.(model[:Rp]).data
    RmHat = value.(model[:Rm]).data

    grp_iter = create_group(grp_loop, string(seq_iter))
    if !isnothing(logfile)
      grp_iter["solver.log"] = read(logfile, String)
      rm(logfile)
      logfile = nothing
    else
      grp_iter["solver.log"] = "NA"
    end

    grp_iter["solve.status"] = Int(termination_status(model))
    grp_iter["solve.status.str"] = string(termination_status(model))
    grp_iter["solve.time"] = solve_time(model)
    grp_iter["solve.iter.barrier"] = barrier_iterations(model)
    grp_iter["solve.iter.simplex"] = simplex_iterations(model)
    grp_iter["solver"]  = solver



    if !isnothing(RpHat_old) && !isnothing(RmHat_old)
      errp = norm((RpHat - RpHat_old)./RpHat_old, norm_err)
      errm = norm((RmHat - RmHat_old)./RmHat_old, norm_err)
      iter_str = rpad("Iter $(seq_iter): ", 12)
      offset_str = rpad("", 12)
      @info BLUE_FG("$(iter_str)linearization error (ϵR+, ϵR-) = ($errp,$errm)")
      grp_iter["Rp.err"] = errp
      grp_iter["Rm.err"] = errm
      p_old, flux_old = Riemann2physics(RpHat_old, RmHat_old, c)
      p, flux = Riemann2physics(RpHat, RmHat, c)
      err_pressure = norm((p - p_old)./p_old, norm_err)
      err_flux = norm((flux - flux_old)./flux_old, norm_err)
      @info BLUE_FG("$(offset_str)linearization error (ϵp, ϵflux) = ($err_pressure,$err_flux)")
      grp_iter["p.err"] = err_pressure
      grp_iter["q.err"] = err_flux

      # Try to see if we cycle
      if !isnothing(RpHat_oldold) && !isnothing(RmHat_oldold)
        errp2 = norm((RpHat - RpHat_oldold)./RpHat_oldold, norm_err)
        errm2 = norm((RmHat - RmHat_oldold)./RmHat_oldold, norm_err)
        @info BLUE_FG("$(iter_str)linearization error (ϵ2R+, ϵ2R-) = ($errp2,$errm2)")
        grp_iter["Rp.err.2"] = errp2
        grp_iter["Rm.err.2"] = errm2
        p_oldold, flux_oldold = Riemann2physics(RpHat_oldold, RmHat_oldold, c)
        p2, flux2 = Riemann2physics(RpHat, RmHat, c)
        err_pressure2 = norm((p2 - p_oldold)./p_oldold, norm_err)
        err_flux2 = norm((flux2 - flux_oldold)./flux_oldold, norm_err)
        @info BLUE_FG("$(offset_str)linearization error (ϵ2p, ϵ2flux) = ($err_pressure2,$err_flux2)")
        grp_iter["p.err.2"] = err_pressure2
        grp_iter["q.err.2"] = err_flux2

      end
      RpHat_oldold = RpHat_old
      RmHat_oldold = RmHat_old

      # break if we are less than the tolerance
      if errp < seq_err_tol && errm < seq_err_tol
        break
      end

    end

    if (seq_iter >= seq_iter_max)
      @info BLUE_FG("Reached max sequential iter $(seq_iter_max)")
      break
    end

    RpHat_old = RpHat
    RmHat_old = RmHat

    params["FricL"] = get_linearized_friction_term_from_Riemann(RpHat, RmHat)


    # Prepare next iteration
    model, logfile = get_model(solver, solvers_opts)
    seq_iter += 1
  end

  # dictionary of numerical values
  vals = Dict{String, Any}()
  vals["objective_value"] = objective_value(model) / num_coeffs["scaling_obj"] 
  vals["penalization_terminal_state_value"] = value(solve_artifacts["terminal_penalization_control_expr"])
  vals["penalization_control_value"] = value(solve_artifacts["penalization_control_expr"])
  vals["profit"] = value(solve_artifacts["profit_expr"])
  vals["seq_iter"] = seq_iter
  vals["seq_err:R"] = (errp, errm)
  vals["seq_err:pq"] = (err_pressure, err_flux)

  attrs = attributes(fid)
  attrs["nb_time_pts"] = nb_t_pts
  attrs["nb_x_pts"] = nb_x_pts
  attrs["T_c"] = T_c
  attrs["save_multipliers"] = save_multipliers
  attrs["delta_t"] = params["delta_t"]
  attrs["delta_x"] = params["delta_x"]

  fid["optimal_stationary_flows"] = solve_artifacts["opt_mflows"]
  fid["penalization_value"] = vals["penalization_control_value"]
  fid["penalization_terminal_state"] = vals["penalization_terminal_state_value"]
  fid["profit"] = vals["profit"]

  attrs["solver"] = solver
  attrs["solvers_opts"] = JSON.json(solvers_opts)
  attrs["model_opts"] = JSON.json(model_opts)

  attrs["num_vars"] = num_variables(model)
  cons_types = list_of_constraint_types(model)
  attrs["num_cons"] =  sum(num_constraints(model, expr, type) for (expr, type) in cons_types)

  attrs["seq_iter"] = seq_iter

  if !isnothing(logfile)
    attrs["logfile"] = read(logfile, String)
    rm(logfile)
  else
    attrs["logfile"] = "NA"
  end

  save_data(model, fid, params, save_multipliers=save_multipliers)

  return (model, hdf5_filename, vals)
end
