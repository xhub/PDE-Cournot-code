# SPDX-License-Identifier: MIT
using ImageMagick, Gurobi

include("hPDE_DO.jl")

if length(ARGS) != 1
  println("ERROR: you need to give only one image file as argument")
end



dat = magickinfo(ARGS[1], "jsons")
jsons = dat["jsons"]

if isnothing(jsons)
  println("ERROR: File $(ARGS[1]) is missing the exif attribute \"jsons\"")
end

for (k, v) in  JSON.parse(jsons)
  println("Reproducing file $k")
  nb_x_pts = v["nb_x_pts"]
  solver = v["solver"]
  solvers_opts = JSON.parse(v["solvers_opts"])
  model_opts = JSON.parse(v["model_opts"])
  cfg_json = JSON.parse(v["cfg_json"])

  json_data = Dict{String,Any}(key => cfg_json[replace(key, r".json$"=>"")]
                               for key in (k * ".json" for k in keys(cfg_json)))

  mm, hdf5_fname, vals = hPDE_DO(json_data=json_data,
            nb_x_pts=nb_x_pts,
            solvers_opts=solvers_opts,
            model_opts=model_opts,
            solver=solver,
            cfg_dir="/dev/null") # this is just to ensure that we don't misread the data

  println("Reproduced file is $hdf5_fname")
end
