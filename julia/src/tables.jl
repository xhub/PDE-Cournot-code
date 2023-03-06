# SPDX-License-Identifier: MIT
using ArgParse, HDF5, JSON, Printf

include("hPDE_DO.jl")

function _key2float(key::String, params)
  red, pmax = params
  if key == "pmax"
    return pmax
  elseif key == "red"
    return red
  end

  error("Unsupported key $(key)")
end

function _store_vals(key, profits, penalization_rel_vals, prod_surplus, cons_surplus, profit, pen_rel_val, ps, cs)
	profits[key] = profit
	penalization_rel_vals[key] = pen_rel_val
	prod_surplus[key] = ps
	cons_surplus[key] = cs
end

function store_vals(key::String, params, profits, penalization_rel_vals, prod_surplus, cons_surplus, profit, pen_rel_val, ps, cs)
  numkey = _key2float(key, params)
  _store_vals(numkey, profits, penalization_rel_vals, prod_surplus, cons_surplus, profit, pen_rel_val, ps, cs)
end

function _r(v::Float64)::String
  return @sprintf "%1.2e" v
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--index"
            help = "value used as row parameter"
            arg_type = String
            required = true
            default = "pmax"
        "files"
            help = "HDF5 files to process"
            nargs = '+'
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

key = parsed_args["index"]

reds = Vector{Float64}()
pmaxs = Vector{Float64}()
profits =  Dict{Float64,Float64}()
penalization_rel_vals = Dict{Float64,Float64}()
prod_surplus = Dict{Float64,Float64}()
cons_surplus = Dict{Float64,Float64}()
idx2file = Dict{Float64,String}()
pmax2file = Dict{Float64,String}()
red2file = Dict{Float64,String}()

for arg in parsed_args["files"]
	params, profit, pen_rel_val, ps, cs = h5open(arg) do f

		attrs = attributes(f)
		cfg_json = JSON.parse(read(attrs["cfg_json"]))
		solvers_opts = JSON.parse(read(attrs["solvers_opts"]))
		model_opts = JSON.parse(read(attrs["model_opts"]))
		json_data = Dict{String,Any}(key => cfg_json[replace(key, r".json$"=>"")]
                               for key in (k * ".json" for k in keys(cfg_json)))
		econ_data = cfg_json["econ"]
		pipe_data = cfg_json["pipe"]

		# Get the upper pressure bound
		pmax = pipe_data["pressure_max"]/1e5

		# Compute the reduction value
		cost_1 = econ_data["cost_functions"][1]["intercept"]
		cost_2 = econ_data["cost_functions"][2]["intercept"]

		red_cost = 1 - cost_1/cost_2

		obj_value = read(f["obj_value"])
		penalization = read(f["penalization_terminal_state"]) + read(f["penalization_value"])
		relative_pen_val = penalization / (obj_value)

		profit = read(f["profit"])

		qin = read(f["control"]["mflowIN"])
		qout = read(f["control"]["mflowOUT"])
		Δt = read(attrs["delta_t"])
		T_c = read(attrs["T_c"])
		times = collect(range(0., T_c; length=read(attrs["nb_time_pts"])))

		@assert abs(times[2] - times[1] - Δt) < 1e-10

		# Compute the values
		cs = consumer_surplus_discrete(qout, times, econ_data)
		ps = producer_surplus(qin, qout, times, econ_data)

		println("\"profit\" - welfare: $(profit - ps - cs)")
    params = (red_cost, pmax)
		return (params, profit, relative_pen_val, ps, cs)
	end

  red, pmax = params

	push!(reds, red)
  push!(pmaxs, pmax)

  store_vals(key, params, profits, penalization_rel_vals, prod_surplus, cons_surplus, profit, pen_rel_val, ps, cs)
	pmax2file[pmax] = arg
	red2file[red] = arg
end

if key == "pmax"
  row_idx = pmaxs
  row_label = "\\bar{p}"
  file_dict = pmax2file
elseif key == "red"
  row_idx = reds
  row_label = "reduction"
  file_dict = red2file
end

println("$(row_label) & profit & penalization & producer surplus & consumer surplus & welfare")
for row in sort(row_idx)
  println("$(row) & $(_r(profits[row])) & $(_r(penalization_rel_vals[row])) & $(_r(prod_surplus[row])) & $(_r(cons_surplus[row])) & $(_r(prod_surplus[row] + cons_surplus[row]))  % $(file_dict[row])")
end
