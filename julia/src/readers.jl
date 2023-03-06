# SPDX-License-Identifier: MIT
function anon_reader(fname::String)
    d = JSON.parsefile(fname)
    return d
end

pipe_reader = anon_reader
econ_reader = anon_reader
approx_reader = anon_reader

