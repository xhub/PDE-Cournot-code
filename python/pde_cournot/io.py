# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

import h5py
import json
import numpy as np
from collections import namedtuple
from collections.abc import Iterable

from .data import Pipe, Econ, Approx


def approx_reader(fname):
    with open(fname, 'r') as f:
        approx_dict = json.load(f)
        return Approx(**approx_dict)


def econ_reader(fname):
    with open(fname, 'r') as f:
        econ_dict = json.load(f)
        return Econ(**econ_dict)


def pipe_reader(fname):
    with open(fname, 'r') as f:
        pipe_dict = json.load(f)
        return Pipe(**pipe_dict)


def json_writer(fname: str, nt: namedtuple):
    with open(fname, 'w') as f:
        json.dump(nt._asdict(), f, indent=2)


pipe_writer = json_writer
econ_writer = json_writer
approx_writer = json_writer


def hdf5_writer_init(fname: str, pipe, econ, approx, solver_cfg=None, version: int=0):
    """ Initialize the HDF5 file for saving the simulation result, and return the file descriptor """
    f = h5py.File(fname, 'w')

    f.attrs["pipe_json"] = json.dumps(pipe._asdict())
    f.attrs["econ_json"] = json.dumps(econ._asdict())
    f.attrs["approx_json"] = json.dumps(approx._asdict())
    if solver_cfg:
        f.attrs["solver_json"] = json.dumps(solver_cfg._asdict())
    f.attrs["version"] = version

    return f

def _hdf5_to_nt(toplevel, attr_name: str, nt: namedtuple):
    if attr_name in toplevel.attrs:
        nt_data = json.loads(toplevel.attrs[attr_name])
        return nt(**nt_data)
    else:
        print(f"Warning: HDF5 has no toplevel attribute name {attr_name}, could not create object {nt}")
        return None


def hdf5_configs_reader(fd: h5py.File):
    """ read configs from an HDF5 file """
    toplevel = fd["/"]
    version = toplevel.attrs['version']
    out = {}
    if version < 3:
        all_data = zip((Pipe, Econ, Approx), ("pipe_json", "econ_json", "approx_json"))
        for (nt, attr_name) in all_data:
            out[attr_name] = _hdf5_to_nt(toplevel, attr_name, nt)
    else:
        cfg_json = json.loads(toplevel.attrs['cfg_json'])
        all_data = zip((Pipe, Econ, Approx), ("pipe_json", "econ_json", "approx_json"))
        for (nt, attr_name) in all_data:
            out[attr_name] = nt(**cfg_json[attr_name.replace('_json', '')])


    return out

def hdf5_save(h5handle, k: str, v):
    """ save a dictionary entry (k,v) to hdf5) """

    type_v = type(v)
    if isinstance(v, float):
        h5handle[k] = v
    elif type_v == dict:
        subfolder = h5handle.create_group(k)
        for kk, vv in v.items():
            hdf5_save(subfolder, kk, vv)
    elif isinstance(v, np.ndarray):
        h5handle[k] = v
    elif isinstance(v, list):
        for (idx, e) in enumerate(v):
            hdf5_save(h5handle, str(idx), e)
    else:
        raise RuntimeError(f"unsupported type {type_v} as value for key {k}")


