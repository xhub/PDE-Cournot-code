# SPDX-License-Identifier: MIT
using HDF5, JSON, Gurobi

include("hPDE_DO.jl")

if length(ARGS) != 1
  println("ERROR: you need to give only one HDF5 file as argument")
  exit()
end

json_data, solvers_opts, model_opts, solver, nb_x_pts = h5open(ARGS[1]) do f
  attrs = attributes(f)
  cfg_json = JSON.parse(read(attrs["cfg_json"]))
  solvers_opts = JSON.parse(read(attrs["solvers_opts"]))
  model_opts = JSON.parse(read(attrs["model_opts"]))
  solver = read(attrs["solver"])
  nb_x_pts = read(attrs["nb_x_pts"])
  json_data = Dict{String,Any}(key => cfg_json[replace(key, r".json$"=>"")]
                               for key in (k * ".json" for k in keys(cfg_json)))
  return (json_data, solvers_opts, model_opts, solver, nb_x_pts)
end

println(keys(json_data))

mm, hdf5_fname, vals = hPDE_DO(json_data=json_data,
            nb_x_pts=nb_x_pts,
            solvers_opts=solvers_opts,
            model_opts=model_opts,
            solver=solver,
            cfg_dir="/dev/null") # this is just to ensure that we don't misread the data


