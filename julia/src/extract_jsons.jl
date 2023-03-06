# SPDX-License-Identifier: MIT
using HDF5, JSON

# macros are magic
macro write_json(json_name)
  name = esc(json_name)
  name_str = string(json_name)
  return :(open("$($name_str).json", "w") do io
             write(io, json($name, 2))
           end)
end

f = h5open(ARGS[1])
attrs = attributes(f)
cfg_json = JSON.parse(read(attrs["cfg_json"]))
solvers_opts = JSON.parse(read(attrs["solvers_opts"]))
model_opts = JSON.parse(read(attrs["model_opts"]))
solver = read(attrs["solver"])
nb_x_pts = read(attrs["nb_x_pts"])

# ease of saving with macro
approx = cfg_json["approx"]
econ = cfg_json["econ"]
pipe = cfg_json["pipe"]

# Tidy things up
# This should no longer be needed on recent runs, but it's good for boostrapping
# and backward compatibility
solvers_opts["current_solver"] = solver
model_opts["nb_x_pts"] = nb_x_pts

@write_json(solvers_opts)
@write_json(model_opts)
@write_json(approx)
@write_json(econ)
@write_json(pipe)
