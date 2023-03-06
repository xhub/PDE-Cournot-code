# SPDX-License-Identifier: MIT
import numpy as np

def _get_data(time, times, objs):
    if len(objs) == 1:
        return objs[0]

    times_len = times
    for i in range(times_len-1):
    # the case with the lsat time interval is special
        if time >= times[i] and (time < times[i+1] or (i+1 == times_len and time == times[i+1])):
            return objs[i]

    raise RuntimeError(f"given time ${time} is outside of the interval [{times[0]}, {times[1]}]")

def get_cost_function(econdata, time):
    """ Return the cost function at a given time"""
    return _get_data(time, econdata.cost_functions, econdata.cost_functions)

def get_inverse_demands(econdata, time):
    """ Return the inverse demand function at a given time """
    return _get_data(time, econdata.time_intervals, econdata.inverse_demands)

def h5attr2comment(v):
    if type(v) is np.bytes_:
        return v.decode('UTF-8')
    elif isinstance(v, np.integer):
        return int(v)
    else:
        return v

def save_model(f) -> dict:
    version = f.attrs['version']
    d = {}
    if version < 3:
        for k, v in f.attrs.items():
            if '_json' in k:
                d[k] = h5attr2comment(v)
    else:
        for k in ('nb_x_pts', 'solver', 'model_opts', 'solvers_opts', 'cfg_json'):
            d[k] = h5attr2comment(f.attrs[k])

    return d

