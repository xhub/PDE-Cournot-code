# SPDX-License-Identifier: MIT
using HDF5, JSON, Gurobi

include("hPDE_DO.jl")

if length(ARGS) < 1
  println("ERROR: you need to give at least one directory as argument")
  exit()
end

function get_complete_cfg(dirname)
   cfg_json = JSON.parsefile(joinpath(dirname, "cfg.json"))
   solvers_opts = JSON.parsefile(joinpath(dirname, "solvers_opts.json"))
   model_opts = JSON.parsefile(joinpath(dirname, "model_opts.json"))
   solver = read(joinpath(dirname, "solver"), String)
   nb_x_pts = parse(Int, read(joinpath(dirname, "nb_x_pts"), String))
   json_data = Dict{String,Any}(key => cfg_json[replace(key, r".json$"=>"")]
                               for key in (k * ".json" for k in keys(cfg_json)))
  return (json_data, solvers_opts, model_opts, solver, nb_x_pts)
end

for arg in ARGS
   json_data, solvers_opts, model_opts, solver, nb_x_pts = get_complete_cfg(arg)
   hPDE_DO(json_data=json_data,
            nb_x_pts=nb_x_pts,
            solvers_opts=solvers_opts,
            model_opts=model_opts,
            solver=solver,
            cfg_dir="/dev/null") # this is just to ensure that we don't misread the data
end

