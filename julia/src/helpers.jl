# SPDX-License-Identifier: MIT
myeps = 1e-6

function get_time_period(econ_data, time::Real)
    if time < econ_data["time_intervals"][1]
        initial_time = econ_data["time_intervals"][1]
        throw(DomainError("time $(time) smaller than initial time $(initial_time)"))
    end
    for (idx, tk) in enumerate(econ_data["time_intervals"][2:end])
        if (tk > time-myeps)
            return idx
        end
    end

    return length(econ_data["time_intervals"])-1
end

function get_obj_csts(econ_data, tk::Real, tkp1::Real)
    period_length = tkp1 - tk
    # WARNING: here we assume that the economical period are multiples of the
    # discretization one
    time_period_idx = get_time_period(econ_data, tkp1)
    beta = econ_data["inverse_demands"][time_period_idx]["intercept"]
    a = econ_data["inverse_demands"][time_period_idx]["slope"]

    if length(econ_data["cost_functions"]) == 1
        cost = econ_data["cost_functions"][1]["intercept"]
    else
        cost = econ_data["cost_functions"][time_period_idx]["intercept"]
    end

    @assert beta > 0 && a < 0

    factor_first_term = period_length * beta
    N = econ_data["number_players"]
    coeff_potential = 1.

    if N > 0
      coeff_potential = (N+1) / (N)
    end
    factor_second_term = period_length * coeff_potential * (a / 2.)


    return (factor_first_term, factor_second_term, cost * period_length)
end

function _getecondata(econ_data, t::Real, smoothing_delta::Real)
    time_period_idx = get_time_period(econ_data, t)
    time_period_ahead  = get_time_period(econ_data, t+smoothing_delta)

    if t-smoothing_delta < 0.
        time_period_behind = 1
    else
        time_period_behind = get_time_period(econ_data, t-smoothing_delta)
    end

    if time_period_ahead > time_period_behind + 1
       @error "time period issue"
    end

    time_period_max = length(econ_data["time_intervals"])

    if time_period_ahead == time_period_behind

        time_period_idx = time_period_ahead
        revenue_cst = econ_data["inverse_demands"][time_period_idx]["intercept"]
        revenue_lin = econ_data["inverse_demands"][time_period_idx]["slope"]
        prodcost_cst = econ_data["cost_functions"][time_period_idx]["intercept"]
        prodcost_lin = econ_data["cost_functions"][time_period_idx]["slope"]

    else
        # In this case, we are close to a switch and we smoothen the data


        revenue_cst_prev = econ_data["inverse_demands"][time_period_behind]["intercept"]
        revenue_cst_succ = econ_data["inverse_demands"][time_period_ahead]["intercept"]
        revenue_lin_prev = econ_data["inverse_demands"][time_period_behind]["slope"]
        revenue_lin_succ = econ_data["inverse_demands"][time_period_ahead]["slope"]

        @assert length(econ_data["cost_functions"]) > 1

        prodcost_cst_prev = econ_data["cost_functions"][time_period_behind]["intercept"]
        prodcost_cst_succ = econ_data["cost_functions"][time_period_ahead]["intercept"]
        prodcost_lin_prev = econ_data["cost_functions"][time_period_behind]["slope"]
        prodcost_lin_succ = econ_data["cost_functions"][time_period_ahead]["slope"]

        tmid = econ_data["time_intervals"][time_period_ahead]
        tlo = tmid - smoothing_delta
        δt = t - tlo
        if δt > 2*smoothing_delta
            @error "δt = $(δt) > $(2*smoothing_delta) = smoothing_delta"
            exit()
        elseif δt < 0
            @error "δt = $(δt)  < 0"
            exit()
        end

        revenue_cst_slope  = (revenue_cst_succ - revenue_cst_prev)/(2*smoothing_delta)
        revenue_lin_slope  = (revenue_lin_succ - revenue_lin_prev)/(2*smoothing_delta)
        prodcost_cst_slope = (prodcost_cst_succ - prodcost_cst_prev)/(2*smoothing_delta)
        prodcost_lin_slope = (prodcost_lin_succ - prodcost_lin_prev)/(2*smoothing_delta)

        revenue_cst  =  revenue_cst_prev  + δt*revenue_cst_slope
        revenue_lin  =  revenue_lin_prev  + δt*revenue_lin_slope
        prodcost_cst =  prodcost_cst_prev + δt*prodcost_cst_slope
        prodcost_lin =  prodcost_lin_prev + δt*prodcost_lin_slope

    end

    return [revenue_cst, revenue_lin, prodcost_cst, prodcost_lin]
end

function get_price(econ_data, tk::Real, qk::Real)::Real
  time_period_idx = get_time_period(econ_data, tk)
  revenue_cst_single_tk = econ_data["inverse_demands"][time_period_idx]["intercept"]
  revenue_lin_single_tk = econ_data["inverse_demands"][time_period_idx]["slope"]

  return revenue_cst_single_tk + revenue_lin_single_tk * qk
end

function get_cost(econ_data, tk::Real, qin::Real)::Real
  time_period_idx = get_time_period(econ_data, tk)
  if length(econ_data["cost_functions"]) == 1
    prodcost_cst_single_tk = econ_data["cost_functions"][1]["intercept"]
    prodcost_lin_single_tk = econ_data["cost_functions"][1]["slope"]
  else
    prodcost_cst_single_tk = econ_data["cost_functions"][time_period_idx]["intercept"]
    prodcost_lin_single_tk = econ_data["cost_functions"][time_period_idx]["slope"]
  end

  return prodcost_cst_single_tk + prodcost_lin_single_tk * qin
end

function get_obj_csts_v2(econ_data, tk::Real, tkp1::Real; smoothing_delta::Real=0.)
    period_length = tkp1 - tk

    # WARNING: here we assume that the economical period are multiples of the
    # discretization one

    if smoothing_delta < 0.0
        @error "smoothing parameter is negative"
        exit()

    elseif smoothing_delta > 0.0

        revenue_cst_single_tk, revenue_lin_single_tk, prodcost_cst_single_tk, prodcost_lin_single_tk = _getecondata(econ_data, tk, smoothing_delta)
        revenue_cst_single_tkp1, revenue_lin_single_tkp1, prodcost_cst_single_tkp1, prodcost_lin_single_tkp1 = _getecondata(econ_data, tkp1, smoothing_delta)

    else
        inv_demands = econ_data["inverse_demands"]
        cost_fns    = econ_data["cost_functions"]

        time_period_idx_k   = get_time_period(econ_data, tk)
        time_period_idx_kp1 = get_time_period(econ_data, tkp1)

        revenue_cst_single_tk   = inv_demands[time_period_idx_k]["intercept"]
        revenue_cst_single_tkp1 = inv_demands[time_period_idx_kp1]["intercept"]
        revenue_lin_single_tk   = inv_demands[time_period_idx_k]["slope"]
        revenue_lin_single_tkp1 = inv_demands[time_period_idx_kp1]["slope"]

        if length(cost_fns) == 1
            prodcost_cst_single_tk = prodcost_cst_single_tkp1 = cost_fns[1]["intercept"]
            prodcost_lin_single_tk = prodcost_lin_single_tkp1 = cost_fns[1]["slope"]
        else
            prodcost_cst_single_tk   = cost_fns[time_period_idx_k]["intercept"]
            prodcost_cst_single_tkp1 = cost_fns[time_period_idx_kp1]["intercept"]
            prodcost_lin_single_tk   = cost_fns[time_period_idx_k]["slope"]
            prodcost_lin_single_tkp1 = cost_fns[time_period_idx_kp1]["slope"]
        end
    end

    # By now, revenue_cst_single_tk, revenue_cst_single_tkp1, prodcost_cst_single_tk,
    #         prodcost_cst_single_tkp1, revenue_lin_single_tk, revenue_lin_single_tkp1,
    #         prodcost_lin_single_tk, prodcost_lin_single_tkp1
    # are well defined

    @assert revenue_cst_single_tk > 0 && revenue_cst_single_tkp1 > 0 && revenue_lin_single_tk < 0 && revenue_lin_single_tkp1 < 0

    revenue_cst_tk    = period_length * revenue_cst_single_tk
    revenue_cst_tkp1  = period_length * revenue_cst_single_tkp1
    prodcost_cst_tk   = period_length * prodcost_cst_single_tk
    prodcost_cst_tkp1 = period_length * prodcost_cst_single_tkp1

    # Small hack to support the welfare case
    N = econ_data["number_players"]
    coeff_potential = 1.
    if N > 0
      coeff_potential = (N+1) / (N)
    end

    #                 discretization * coeff from potential * 1/2 [from quadratic]
    factor_lin_term   = period_length   * coeff_potential   * 1 / 2

    revenue_lin_tk    = factor_lin_term * revenue_lin_single_tk
    revenue_lin_tkp1  = factor_lin_term * revenue_lin_single_tkp1
    prodcost_lin_tk   = factor_lin_term * prodcost_lin_single_tk
    prodcost_lin_tkp1 = factor_lin_term * prodcost_lin_single_tkp1

    econ_csts = Dict{String,Float64}(
      "revenue_cst_tk"    => revenue_cst_tk,
      "revenue_cst_tkp1"  => revenue_cst_tkp1,
      "revenue_lin_tk"    => revenue_lin_tk,
      "revenue_lin_tkp1"  => revenue_lin_tkp1,
      "prodcost_cst_tk"   => prodcost_cst_tk,
      "prodcost_cst_tkp1" => prodcost_cst_tkp1,
      "prodcost_lin_tk"   => prodcost_lin_tk,
      "prodcost_lin_tkp1" => prodcost_lin_tkp1,
    )


    return econ_csts
end

function physic2Riemann(p, flux, c)
    M = c * flux ./ p
    Rp = (M .+ 1) .* p
    Rm = (M .- 1) .* p
    return (Rp, Rm)
end

function Riemann2physics(Rp, Rm, c)
    p = (Rp - Rm) / 2
    flux = (Rp + Rm) / (2 * c)
    return (p,flux)
end

function p_init_LISO2(pOP::Float64, qOP::Float64, c::Float64, xs::Vector{Float64}, pipe_data::Dict)
    p0 = pOP
    D = pipe_data["diameter"]
    λ = nikuradse(D, pipe_data["k"])
    σ = sigma(D, λ, c)

    α = σ * qOP^2 / (2*pOP^2)

    return p0 * (2 .- exp.(α * xs))
end

# Basic sanity check: is the initial state feasible?

function check_initial_state_feasible_LISO2(pipe_data, approx_data)
    D = pipe_data["diameter"]
    λ = nikuradse(D, pipe_data["k"])
    speed_of_sound = sqrt(pipe_data["T"] * pipe_data["Rs"])
    L = pipe_data["length"]
    σ = 16 * speed_of_sound^2 * λ / (π^2 * D^5)
    qOP = approx_data["operational_massflow"]
    p0 = pOP = approx_data["operational_pressure"]

    α = σ * qOP^2 / (2*pOP^2)
    pL = p0 * (2 - exp(α * L))

    @info "LISO2: pin = $(p0); pout = $(pL); mflow = $(qOP)"

    if pL <=  pipe_data["pressure_min"]
      @error "pout = $(pL) < $(pipe_data["pressure_min"]) the lower bound"
      exit()
    end
end

function check_initial_state_feasible_ISO2(pipe_data, approx_data)
    friction = nikuradse(pipe_data["diameter"], pipe_data["k"])
    speed_of_sound = sqrt(pipe_data["T"] * pipe_data["Rs"])
    Λ = (16 * speed_of_sound^2 * friction)/(π^2*pipe_data["diameter"]^5) * pipe_data["length"]
    mflow = approx_data["operational_massflow"]
    pin  = approx_data["operational_pressure"]
    pout = sqrt(pin^2 - Λ*mflow^2)

    @info "ISO2: pin = $(pin); pout = $(pout); mflow = $(mflow); Λ = $(Λ)"

    if pout <=  pipe_data["pressure_min"]
      @error "pout = $(pout) < $(pipe_data["pressure_min"]) the lower bound"
      exit()
    end
end

function get_num_coeffs_default(Δx::Float64, Δt::Float64)
    # The default value give a coefficient around 1 for the unknowns
    num_coeffs = Dict{String,Any}("penalty_terminal_curvature"  =>  Δx^4,
                                  "penalty_terminal_gradient"   =>  Δx^2,
                                  "penalty_obj_ctrl_laplacian"  =>  Δt^4/Δt,
                                  "penalty_obj_ctrl_grad"       =>  Δt^2,
                                  "penalty_obj_ctrl_grad_p"     =>  Δt^2/Δt * 1e-7,
                                  "penalty_obj_ctrl_grad_q"     =>  Δt^2/Δt * 1e-2,
                                  "penalty_obj_ctrl_L2"         =>  0.,
                                  "scaling_obj"                 =>  1e-1,
                                 )

    return num_coeffs
end

function get_num_coeffs_default(pipe_data, nb_x_pts::Int)
    Δx = pipe_data["length"] / (nb_x_pts - 1)
    c = sqrt(pipe_data["T"] * pipe_data["Rs"])
    Δt = Δx * c

    return get_num_coeffs_default(Δx, Δt)
end

function get_linearized_friction_term(pOP::Float64, fluxOP::Float64, c::Float64, scale_R::Float64)
# We impose the flow direction
  RpOP, RmOP = physic2Riemann(pOP, fluxOP, c)
  @assert RpOP - RmOP > 0.
  @assert RpOP + RmOP > 0.

  # This is the affine approximation of
  # F(Rp, Rm) = (Rp + Rm)^2 / (Rp - Rm)
  # As used in the paper with Gugat, Habermann, Hintermüller and Huber
  # Here we assume the flow to be positive to get rid of the absolute value.
  #
  nabla_F_Rp = (RpOP - 3*RmOP)*(RpOP + RmOP)/(RpOP - RmOP)^2
  nabla_F_Rm = (3*RpOP - RmOP)*(RpOP + RmOP)/(RpOP - RmOP)^2
  FvalOP = (RpOP + RmOP)^2/(RpOP - RmOP)

  # we need to scale the operational point
  FricL(m,Rp,Rm, tidx, xidx) = @expression(m, nabla_F_Rp * (Rp[tidx, xidx] - scale_R * RpOP)
                                            + nabla_F_Rm * (Rm[tidx, xidx] - scale_R * RmOP)
                                            + scale_R * FvalOP)
end

function get_linearized_friction_term_from_Riemann(RpHat::Matrix{Float64}, RmHat::Matrix{Float64})
# We impose the flow direction
#  RpHat - RmHat > 0.
#  RpHat + RmHat > 0.

  # This is the affine approximation of
  # F(Rp, Rm) = (Rp + Rm)^2 / (Rp - Rm)
  # As used in the paper with Gugat, Habermann, Hintermüller and Huber
  # Here we assume the flow to be positive to get rid of the absolute value.
  #
  nabla_F_Rp = (RpHat .- 3*RmHat).*(RpHat .+ RmHat)./(RpHat .- RmHat).^2
  nabla_F_Rm = (3*RpHat .- RmHat).*(RpHat .+ RmHat)./(RpHat .- RmHat).^2
  FvalHat = ((RpHat .+ RmHat).^2) ./(RpHat .- RmHat)

  # we need to scale the operational point
  FricL(m, Rp, Rm, tidx, xidx) = @expression(m, nabla_F_Rp[1+tidx, 1+xidx] * (Rp[tidx, xidx] - RpHat[1+tidx, 1+xidx])
                                              + nabla_F_Rm[1+tidx, 1+xidx] * (Rm[tidx, xidx] - RmHat[1+tidx, 1+xidx])
                                              + FvalHat[1+tidx, 1+xidx])
end
