#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

from collections import OrderedDict
from colorama import init, Fore, Back, Style
import copy
from datetime import datetime
import h5py
from itertools import cycle
import json
import matplotlib
matplotlib.use("Agg")
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.animation as manimation
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import pathlib
import re
import sys

# local imports
import pde_cournot
pde_cournot.set_latex_preamble()

# colorama
init()

#LS = cycle(["solid", "dashed", "dashdot", "dotted",
#            (0, (3, 1, 1, 1)),       # densely dashdotted
#            (0, (5, 1)),             # densely dashed
#            (0, (1, 1)),             # densely dotted
#            (0, (3, 1, 1, 1, 1, 1)), # densely dash dotdotted
#           ])

LS = cycle([
            (0, (3, 1, 1, 1, 1, 1)), # densely dash dotdotted
            (0, (1, 1)),             # densely dotted
            (0, (5, 1)),             # densely dashed
            (0, (3, 1, 1, 1)),       # densely dashdotted
            "dotted",
           ])

color = cycle([None])

# was (16,7.5)
figsize = (12,4.5)
# width/height
figsize_single = (6,4.5)


def printred(s: str):
    print(Fore.RED + s + Style.RESET_ALL)


def compute_no_network_sol(N: float, cost_fn: dict[str, float], rev_fn: dict[str, float]) -> float:
    if abs(cost_fn["slope"]) > 0:  # lazy test
        print("ERROR: non-zero cost slope")

    cost_cst = cost_fn["intercept"]
    rev_slope = rev_fn["slope"]
    rev_cst = rev_fn["intercept"]

    coeff_potential = 1.
    if N > 0:
        coeff_potential = (N+1)/N

    q = -(rev_cst - cost_cst)/(rev_slope*coeff_potential)
    return q


def save_(d: dict, k: str, v):
    d[k] = pde_cournot.h5attr2comment(v)


def main(hdf5_fnames: list[str], labels: list[str], args):
    #######################################################
    # We follow the following logic:                      #
    # - first we collect all the data from the hdf5 files #
    # - then we create all the figures                    #
    #######################################################
    qINs = {}
    qOUTs = {}
    pINs = {}
    pOUTs = {}
    cfgs = {}
    delta_x = {}
    fnames2label = {}
    nb_t_pts = {}
    pmaxs = {}

    jsons_dict = {}
    endtime: float = float("nan")
    for fname in hdf5_fnames:
        with h5py.File(fname, 'r') as f:
            res = dict(f)

            nb_t_pts[fname] = f.attrs['nb_time_pts']
            delta_x[fname] = f.attrs['delta_x']
            qINs[fname] = np.asarray(res['control']['mflowIN'])
            qOUTs[fname] = np.asarray(res['control']['mflowOUT'])

            pINs[fname] = np.asarray(res['control']['pIN'])
            pOUTs[fname] = np.asarray(res['control']['pOUT'])

            cfg = pde_cournot.hdf5_configs_reader(f)
            cfgs[fname] = cfg
            pmaxs[fname] = cfg['pipe_json'].pressure_max/1e5
#            econ_jsons[fname] = cfg['econ']

            jsons_dict[fname] = pde_cournot.save_model(f)

            if np.isnan(endtime):
                if f.attrs['version'] >= 2:
                    endtime = f.attrs['T_c']
                else:
                    endtime = cfg['econ_json'].time_intervals[-1]

            fnames2label[fname] = pde_cournot.get_label(labels, cfg, args)
#        fnames2label[fname] = r'$\bar{{p}}$ = {:.1f} bar'.format(pipe_jsons[fname]['pressure_max']/1e5)

    label2fnames = {}
    for k, v in fnames2label.items():
        label2fnames[v] = k

    if endtime > 3600 * 2:
        endtime /= 3600
        time_unit = 'h'
    elif endtime > 60*5:
        endtime /= 60
        time_unit = "m"
    else:
        time_unit = "s"

    d = datetime.now()
    date_suffix = d.isoformat(timespec='seconds').replace(':', '_')

    metadata = {'jsons': json.dumps(jsons_dict),
                'label2fnames': json.dumps(label2fnames),
                'cmd line': ' '.join(sys.argv)}

    #############################################################################################
    # Plot mass flow over time
    #############################################################################################

    fig, (axIN, axOUT) = plt.subplots(1, 2, figsize=figsize)  # was (16,7.5)

    axIN.set_title(r"Mass flow IN $q^{\mathrm{in}}$ (kg/s)")
    axOUT.set_title(r"Mass flow OUT $q^{\mathrm{out}}$ (kg/s)")

    axIN.set_xlabel(f"Time ({time_unit})")
    axOUT.set_xlabel(f"Time ({time_unit})")

    econ_data = cfgs[hdf5_fnames[0]]['econ_json']
    N = econ_data.number_players
    opts = [compute_no_network_sol(N, econ_data.cost_functions[i], econ_data.inverse_demands[i]) for i in
            range(len(econ_data.cost_functions))]

    tvec = []
    tvecs = []
    pts = []
    nb_markers = 10

    if args.colors is not None:
        line_colors = cycle(args.colors)
    else:
        line_colors = color

    for (idx, opt) in enumerate(opts):
        if idx+1 >= len(econ_data.time_intervals):
            break
        start_time, end_time = econ_data.time_intervals[idx:idx+2]
        if time_unit == "m":
            start_time /= 60
            end_time /= 60
        elif time_unit == "h":
            start_time /= 3600
            end_time /= 3600

        tvec_last = np.linspace(start_time, end_time, nb_markers)
        tvecs.append(tvec_last)
        tvec.extend(tvec_last.tolist())
        pts.extend([opt]*nb_markers)

    leg_stat_opt = OrderedDict()
    if not args.no_label_no_network:
        p1, = axOUT.plot(tvec, pts,  ls=' ', markersize=6, marker='x', c='k')
        leg_stat_opt['No network'] = p1

    if not args.no_label_no_network and args.label_reduced_cost:
        opt_cross = compute_no_network_sol(N, econ_data.cost_functions[0], econ_data.inverse_demands[1])
        p2, = axOUT.plot(tvec_last, [opt_cross]*nb_markers, ls=' ', markersize=6, marker='o', mfc='none', c='k')
        leg_stat_opt['Reduced cost'] = p2

    assert args.weymouth_sol is None or len(args.weymouth_sol) == len(tvecs)
    if args.weymouth_sol is None:
        args.weymouth_sol = []

    for (idx, sol) in enumerate(args.weymouth_sol):
        p, = axOUT.plot(tvecs[idx], [sol]*nb_markers, ls=' ', markersize=6, marker='1', c='b')
        leg_stat_opt['Weymouth'] = p

    if leg_stat_opt:
        leg_extra = axOUT.legend(leg_stat_opt.values(), leg_stat_opt.keys(), numpoints=1, loc=4)

    lLS = copy.deepcopy(LS)
    lcolor = copy.deepcopy(line_colors)
    for fname in hdf5_fnames:
        tvec = np.linspace(0, endtime, nb_t_pts[fname])
        ls = next(lLS)
        c = next(lcolor)
        axIN.plot(tvec, qINs[fname], label=fnames2label[fname], linestyle=ls, c=c)
        axOUT.plot(tvec, qOUTs[fname], label=fnames2label[fname], linestyle=ls, c=c)

    axIN.legend()
    axOUT.legend()

    if leg_stat_opt:
        axOUT.add_artist(leg_extra)

    fig.savefig(f"q-endpoints_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600,
                metadata=metadata)

    #############################################################################################
    # Plot pressure over time
    #############################################################################################

    fig, (axIN, axOUT) = plt.subplots(1, 2, figsize=figsize)  # was (16,7.5)

    axIN.set_title(r"Pressure IN $p^{\mathrm{in}}$ (bar)")
    axOUT.set_title(r"Pressure OUT $p^{\mathrm{out}}$ (bar)")

    axIN.set_xlabel(f"Time ({time_unit})")
    axOUT.set_xlabel(f"Time ({time_unit})")

    lLS = copy.deepcopy(LS)
    lcolor = copy.deepcopy(line_colors)
    for fname in hdf5_fnames:
        tvec = np.linspace(0, endtime, nb_t_pts[fname])
        ls = next(lLS)
        c = next(lcolor)
        axIN.plot(tvec, pINs[fname], label=fnames2label[fname], linestyle=ls, c=c)
        axOUT.plot(tvec, pOUTs[fname], label=fnames2label[fname], linestyle=ls, c=c)
        axIN.axhline(pmaxs[fname], linestyle='-', color='grey', alpha=0.3)

    axIN.legend()
    axOUT.legend()

    fig.savefig(f"p-endpoints_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600,
                metadata=metadata)

    #############################################################################################
    # Plot linepack over time
    #############################################################################################

    fig, ax = plt.subplots(1, 1, figsize=figsize_single)

    ax.set_title("Evolution of linepack")

    ax.set_xlabel(f"Time ({time_unit})")
    if matplotlib.rcParams['text.usetex']:
        ax.set_ylabel("Variation (in \%)")
    else:
        ax.set_ylabel("Variation (in %)")

    lLS = copy.deepcopy(LS)
    lcolor = copy.deepcopy(line_colors)
    for fname in hdf5_fnames:
        linepack = np.empty((nb_t_pts[fname], ))
        pipe = cfgs[fname]['pipe_json']
        D = pipe.diameter
        T = pipe.T
        Rs = pipe.Rs
        cs = (T*Rs)**.5
        factor = math.pi * D**2 * delta_x[fname] / (4 * cs**2)
        with h5py.File(fname, 'r') as f:
            for i in range(nb_t_pts[fname]):
                i_str = str(i)
                # pressure is stored in bars
                p = 1e5 * np.asarray(f['pipe'][i_str]['pressure'])
                linepack[i] = factor * (np.sum(p[1:-2]) + (p[0] + p[-1])/2)

        initial_linepack = linepack[0]
        linepack_var_perc = 100*(linepack[:]/initial_linepack - 1.)

        tvec = np.linspace(0, endtime, nb_t_pts[fname])
        ls = next(lLS)
        c = next(lcolor)
        np.savetxt('/tmp/tvec.txt', tvec)
        np.savetxt(f'/tmp/linepack_{fname}.txt', linepack_var_perc)
        ax.plot(tvec, linepack_var_perc, label=fnames2label[fname], linestyle=ls, c=c)

    ax.legend()
    fig.savefig(f"linepack-{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600,
                metadata=metadata)


if __name__ == "__main__":
    parser = pde_cournot.get_common_parser()
    parser.add_argument('--no-label-no-network', action='store_true')
    parser.add_argument('--label-reduced-cost', action='store_true')
    parser.add_argument('--weymouth-sol', type=lambda s: [float(item) for item in s.split(',')])

    args = parser.parse_args()

    if not args.label:
        labels = ['pmax']
    else:
        labels = args.label

    hdf5_fnames = args.files

    if not hdf5_fnames:
        printred("No HDF5 file names given!")
        sys.exit(1)

    main(hdf5_fnames, labels, args)
