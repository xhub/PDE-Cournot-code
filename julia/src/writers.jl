# SPDX-License-Identifier: MIT
using HDF5

# Initialize the HDF5 file for saving the simulation result, and return the file descriptor
function hdf5_writer_init(fname, pipe, econ, approx, solver_cfg=nothing)::HDF5.File
    f = h5open(fname, "w")
    attrs = attributes(f)

    attrs["cfg_json"] = JSON.json(Dict{String,Any}("pipe" => pipe, "econ" => econ, "approx" => approx))
    if solver_cfg != nothing
        attrs["solver_json"] = JSON.json(solver_cfg)
    end

    attrs["version"] = 3

    return f
end

function hdf5_init_loop_log(grp::HDF5.H5DataStore, grp_name::String; kwargs...)::HDF5.Group
    grp_loop = create_group(grp, grp_name)
    attrs = attributes(grp_loop)

    for (k, v) in kwargs
        attrs[string(k)] = v
    end

    return grp_loop
end
