# SPDX-License-Identifier: MIT
from collections import namedtuple

pipe_data = ["length", "diameter", "T", "Rs", "slope", "k", "pressure_min", "pressure_max", "flow_min", "flow_max"]
econ_data = ["time_intervals", "number_players", "inverse_demands", "cost_functions"]
approx_data = ["epsilon", "operational_pressure", "operational_massflow"]

Pipe = namedtuple('Pipe', pipe_data)
Econ = namedtuple('Econ', econ_data)
Approx = namedtuple('Approx', approx_data)
