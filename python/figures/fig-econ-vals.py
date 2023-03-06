#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

import argparse
import copy
from datetime import datetime
import glob
import h5py
from itertools import cycle, islice
import json
import matplotlib.pyplot as plt
import numpy as np
import os.path
import re
import sys

# local import
import pde_cournot
pde_cournot.set_latex_preamble()

myeps = 1e-6

#LS = cycle(["solid", "dashed", "dashdot", "dotted",
#            (0, (3, 1, 1, 1)),        # densely dashdotted
#            (0, (5, 1)),              # densely dashed
#            (0, (1, 1)),              # densely dotted
#            (0, (3, 1, 1, 1, 1, 1)),  # densely dash dotdotted
#           ])

LS = cycle([
            (0, (3, 1, 1, 1, 1, 1)),  # densely dash dotdotted
            (0, (1, 1)),              # densely dotted
            (0, (5, 1)),              # densely dashed
            (0, (3, 1, 1, 1)),        # densely dashdotted
            "dotted",
            "dashdot",
            "dashed",
            "solid",
           ]
           )

color = cycle([None])

# was (16,7.5)
figsize = (12, 4.5)
# width/height
figsize_single = (6, 4.5)
lw_surplus = 3

color_PS = 'g'
color_CS = 'b'
color_W = 'k'
color_MPR = 'm'
color_LCS = 'c'
color_NR = 'orange'
color_LRL = 'r'
color_CR = 'lime'

color_price = 'darkred'
color_cost = 'grey'
color_markup = 'navy'

ls_W = "solid"
ls_3 = "dashdot"
ls_10 = "dashed"
ls_20 = "dotted"
ls_30 = (0, (3, 1, 1, 1))


with_stack = False


def get_xaxis_data(labels: list[str], xdicts: dict):
    lab = labels[0]

    if labels[0] == 'pmax':
        xlabel = r"$\bar{p}$ (bar)"
    elif labels[0] == 'qlo':
        xlabel = r"$\underaccent{\bar}{q}$"
    else:
        print("ERROR: unsupported labels value")
        sys.exit(1)

    xvals = xdicts[lab]
    return (xvals, xlabel)


def get_time_period(econ_data: pde_cournot.Econ, time: float):
    if time < econ_data.time_intervals[0]:
        initial_time = econ_data.time_intervals[0]
        raise ValueError(f"time ${time} smaller than initial time ${initial_time}")

    # Here the enumerate is over time intervals with the first one dropped:
    # The returned value is then directly the index
    for (idx, tk) in enumerate(econ_data.time_intervals[1:]):
        if (tk > time-myeps):
            return idx

    return len(econ_data.time_intervals)-2


def get_price(econ_data: pde_cournot.Econ, tk: float, qout: float):
    time_period_idx = get_time_period(econ_data, tk)
    revenue_cst_single_tk = econ_data.inverse_demands[time_period_idx]["intercept"]
    revenue_lin_single_tk = econ_data.inverse_demands[time_period_idx]["slope"]

    return revenue_cst_single_tk + revenue_lin_single_tk * qout


def get_cost(econ_data: pde_cournot.Econ, tk: float, qin: float):
    time_period_idx = get_time_period(econ_data, tk)
    if len(econ_data.cost_functions) == 1:
        prodcost_cst_single_tk = econ_data.cost_functions[1]["intercept"]
        prodcost_lin_single_tk = econ_data.cost_functions[1]["slope"]
    else:
        prodcost_cst_single_tk = econ_data.cost_functions[time_period_idx]["intercept"]
        prodcost_lin_single_tk = econ_data.cost_functions[time_period_idx]["slope"]

    return prodcost_cst_single_tk + prodcost_lin_single_tk * qin


# TODO: does this makes any sense?
# We know that the optimal solution obeys the relation
#
# 0 = \tilde{P}(t, qout(t)) - cost(t) + \tilde{P}'(t, qout(t)) qout(t)
#
# Hence, we can reverse-engineer the cost on each period in the following way:
#
# cost(T_i) = \tilde{P}(T_i, qout(T_i)) + \tilde{P}'(T_i, qout(T_i)) qout(T_i)
#

def compute_mc_nocongestion(qout: list[float], times: np.ndarray, econ_data: pde_cournot.Econ):
    # This is a special task: we are going to infer the marginal cost
    # from the absence of congestion rent.

    time_period_idxs = np.asarray([get_time_period(econ_data, tk) for tk in times])
    time_periods = np.unique(time_period_idxs)

    # Now for each period, we get all the corresponding tks.
    # Then, we retain the middle half of these

    marginal_cost = np.full(times.shape, np.nan)

    # Just for debugging purposes

    for v in time_periods:
        subtkidx = time_period_idxs == v
        subtk = times[subtkidx]
        subqout = qout[subtkidx]
        lq = int(len(subtk)/4)
        tk4avg = subtk[lq:-lq]
        qout4avg = subqout[lq:-lq]

        qoutavg = np.mean(qout4avg)

        Pslope = get_P_slope(econ_data, tk4avg[0])

        num_players = econ_data.number_players
        if num_players == 0:
            inv_num_players = 0
        else:
            inv_num_players = 1./num_players

        marginal_cost[subtkidx] = get_price(econ_data, tk4avg[0], qoutavg) + Pslope*qoutavg*inv_num_players

    return marginal_cost


def get_P_slope(econ_data: pde_cournot.Econ, tk: float) -> float:
    time_period_idx = get_time_period(econ_data, tk)
    revenue_lin_single_tk = econ_data.inverse_demands[time_period_idx]["slope"]

    return revenue_lin_single_tk


def compute_integral(vals: list[tuple], num_expr) -> float:
    Δt = vals[1][0] - vals[0][0]

    return Δt * (sum(num_expr(*v) for v in vals[1:-2]) + num_expr(*vals[0])/2 + num_expr(*vals[-1])/2)


def consumer_surplus(qout: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ qout(t) * [max_price - price(t)] dt

    def num_expr(tk, qk): return (get_price(econ_data, tk, 0.) - get_price(econ_data, tk, qk)) * qk/2
    vals = list(zip(times, qout))

    return compute_integral(vals, num_expr)


def producer_surplus(qin: list[float], qout: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ qout(t) * price(t) - qin(t) * cost(t)] dt
    def num_expr(tk, qink, qoutk): return get_price(econ_data, tk, qoutk)*qoutk - get_cost(econ_data, tk, qink)*qink
    vals = list(zip(times, qin, qout))

    return compute_integral(vals, num_expr)


def revenue_no_network(qout: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ qout(t) * [price(t) - cost(t)] dt

    def num_expr(tk, qoutk): return qoutk * (get_price(econ_data, tk, qoutk) - get_cost(econ_data, tk, qoutk))
    vals = list(zip(times, qout))

    return compute_integral(vals, num_expr)


def total_revenue(qout: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ qout(t) * price(t) dt

    def num_expr(tk, qoutk): return qoutk * get_price(econ_data, tk, qoutk)
    vals = list(zip(times, qout))

    return compute_integral(vals, num_expr)


def linepack_savings(qin: list[float], qout: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ [qin(t) - qout(t)] * cost(t)] dt

    def num_expr(tk, qintk, qoutk): return get_cost(econ_data, tk, qoutk) * (qoutk - qintk)
    vals = list(zip(times, qin, qout))

    return compute_integral(vals, num_expr)


def linepack_loss(qout: list[float], mc: list[float], times: list[float], econ_data: pde_cournot.Econ):
    # ∫ qout(t) * [ mc(t) - cost(t) ] dt

    # TODO: should this be identical?
    #def num_expr(tk, qoutk): return qoutk * (get_price(econ_data, tk, qoutk) - get_cost(econ_data, tk, qoutk))
    def num_expr(tk, qoutk, mctk): return qoutk * (mctk - get_cost(econ_data, tk, qoutk))
    vals = list(zip(times, qout, mc))

    return compute_integral(vals, num_expr)


def congestion_rent(qout: list[float], mc: list[float], times: list[float], econ_data: pde_cournot.Econ) -> float:
    # ∫ qout(t) * [ price(t) - (mc(t) - P_prime * qout(t) / M) ] dt

    num_players = econ_data.number_players
    if num_players == 0:
        inv_num_players = 0
    else:
        inv_num_players = 1./num_players

    def num_expr(tk, qoutk, mctk): return qoutk * ( get_price(econ_data, tk, qoutk) - \
            (mctk - get_P_slope(econ_data, tk) * qoutk * inv_num_players ))
    vals = list(zip(times, qout, mc))

    return compute_integral(vals, num_expr)


def market_power(qout: list[float], times: list[float], econ_data: pde_cournot.Econ) -> float:
    # ∫ qout(t) * [ - P_prime * qout(t) / M ] dt

    num_players = econ_data.number_players
    if num_players == 0:
        inv_num_players = 0
    else:
        inv_num_players = 1./num_players

    def num_expr(tk, qoutk): return -qoutk**2 * get_P_slope(econ_data, tk) * inv_num_players
    vals = list(zip(times, qout))

    return compute_integral(vals, num_expr)


def split_pos_neg(vals):
    return (np.maximum(np.zeros(vals.shape), vals), np.minimum(np.zeros(vals.shape), vals))


def plot_stack_with_neg(ax, xdata: dict, ydata: list[dict]):
    # We do a stack graph allowing for one negative component
    #
    # Credits to
    # https://stackoverflow.com/questions/65859200/how-to-display-negative-values-in-matplotlibs-stackplot

    vals = (d["vals"] for d in ydata)
    labels = list(d["label"] for d in ydata)
    colors = list(d["color"] for d in ydata)

    xvals = xdata["vals"]

    cur_vals = np.zeros(xvals.shape)
    pos_parts = []
    neg_parts = []

    for val in vals:
        cur_vals += val
        pos_part, neg_part = split_pos_neg(cur_vals)
        pos_parts.append(pos_part)
        neg_parts.append(neg_part)

    pos_plot_args = [xvals] + pos_parts
    ax.stackplot(*pos_plot_args, colors=colors, labels=labels)
    # Need to reverse to get things right
    neg_parts.reverse()
    colors.reverse()
    neg_plot_args = [xvals] + neg_parts
    ax.stackplot(*neg_plot_args, colors=colors)

    ax.legend()
    ax.set_xlabel(xdata["label"])


def plot_stack_with_total(ax, xdata: dict, ydata_pos: list[dict], ydata_neg: list[dict],
                          kw_total: dict = {'label': "Firms' profits", 'color': 'g'}):
    # We do a stack graph with positive and negative components
    #
    # Credits to
    # https://stackoverflow.com/questions/65859200/how-to-display-negative-values-in-matplotlibs-stackplot

    vals_pos = list(d["vals"] for d in ydata_pos)
    labels_pos = list(d["label"] for d in ydata_pos)
    colors_pos = list(d["color"] for d in ydata_pos)

    vals_neg = list(d["vals"] for d in ydata_neg)
    labels_neg = list(d["label"] for d in ydata_neg)
    colors_neg = list(d["color"] for d in ydata_neg)

    xvals = xdata["vals"]

    total = np.zeros(xvals.shape)

    for val in vals_pos:
        total += val

    for val in vals_neg:
        total += val

    pos_plot_args = [xvals] + vals_pos
    ax.stackplot(*pos_plot_args, colors=colors_pos, labels=labels_pos)

    if len(vals_neg) > 0:
        neg_plot_args = [xvals] + vals_neg
        ax.stackplot(*neg_plot_args, colors=colors_neg, labels=labels_neg)

    ax.plot(xvals, total, lw=3, **kw_total)

    ax.legend()
    ax.set_xlabel(xdata["label"])


def plot_stack_with_total_and_special(ax, xdata: dict, ydata_pos: list[dict], ydata_special: dict,
                                      kw_total: dict = {'label': "Firms' profits", 'color': 'g'}):
    # We do a stack graph with positive and negative components
    #
    # Credits to
    # https://stackoverflow.com/questions/65859200/how-to-display-negative-values-in-matplotlibs-stackplot

    vals_pos = list(d["vals"] for d in ydata_pos)
    labels_pos = list(d["label"] for d in ydata_pos)
    colors_pos = list(d["color"] for d in ydata_pos)

    vals_special = ydata_special["vals"]
    label_special = ydata_special["label"]
    color_special = ydata_special["color"]

    xvals = xdata["vals"]

    total = np.asarray(copy.deepcopy(vals_special))

    for val in vals_pos:
        total += val

    vals_spos, vals_sneg = split_pos_neg(vals_special)

    pos_plot_args = [xvals, vals_spos]
    pos_plot_args.extend(vals_pos)
    ax.stackplot(*pos_plot_args, colors=[color_special] + colors_pos, labels=[label_special] + labels_pos)

    neg_plot_args = [xvals, vals_sneg]
    ax.stackplot(*neg_plot_args, colors=[color_special])

    ax.plot(xvals, total, lw=3, **kw_total)

    ax.legend()
    ax.set_xlabel(xdata["label"])


def plot_lines_with_total(ax, xdata: dict, ydata_pos: list[dict], ydata_neg: list[dict],
                          kw_plots: dict = {}, kw_total: dict = {'label': "Firms' profits", 'color': 'g'},
                          legend=True):
    # We do a stack graph with positive and negative components
    #
    # Credits to
    # https://stackoverflow.com/questions/65859200/how-to-display-negative-values-in-matplotlibs-stackplot

    vals_pos = list(d["vals"] for d in ydata_pos)
    labels_pos = list(d["label"] for d in ydata_pos)
    colors_pos = list(d["color"] for d in ydata_pos)

    vals_neg = list(d["vals"] for d in ydata_neg)
    labels_neg = list(d["label"] for d in ydata_neg)
    colors_neg = list(d["color"] for d in ydata_neg)

    xvals = xdata["vals"]

    total = np.zeros(xvals.shape)

    for val in vals_pos:
        total += val

    for val in vals_neg:
        total += val

    for vals, c, l in zip(vals_pos, colors_pos, labels_pos):
        ax.plot(xvals, vals, color=c, label=l, **kw_plots)

    for vals, c, l in zip(vals_neg, colors_neg, labels_neg):
        ax.plot(xvals, vals, color=c, label=l, **kw_plots)

    ax.plot(xvals, total, lw=lw_surplus, **kw_total)

    if legend:
        ax.legend()

    ax.set_xlabel(xdata["label"])


def main(fnames: list[str], args, no_plot=False):

    if not args.label:
        labels = ['pmax']
    else:
        labels = args.label

    jsons_dict = {}
    label2fnames = {}
    pmax = []
    qlos = []
    obj = []
    PS = []
    CS = []
    welfare = []

    LCS = []  # linepack savings
    LRL = []  # linepack revenue loss
    MPR = []  # market power rent
    CR = []   # Congestion rent

    qouts = []
    markups = []
    costs = []
    prices = []
    tks = []
    fig_labels = []

    if args.colors is not None:
        line_colors = cycle(args.colors)
    else:
        line_colors = color

    for fname in fnames:
        with h5py.File(fname) as f:
            cfg = pde_cournot.hdf5_configs_reader(f)
            jsons_dict[fname] = pde_cournot.save_model(f)

            econ_data = cfg['econ_json']
            qin = f['control']['mflowIN'][()]
            qout = f['control']['mflowOUT'][()]
            qouts.append(qout)
            T_c = f.attrs["T_c"]
            times = np.linspace(0., T_c, num=f.attrs["nb_time_pts"])
            mc = compute_mc_nocongestion(qout, times, econ_data)

            p_up_bar = cfg['pipe_json'].pressure_max/1e5
            pmax.append(p_up_bar)
            qlo = cfg['pipe_json'].flow_min
            qlos.append(qlo)

            obj.append(f["obj_value"][()])

            PS.append(producer_surplus(qin, qout, times, econ_data))
            CS.append(consumer_surplus(qout, times, econ_data))
            welfare.append(PS[-1] + CS[-1])

            LCS.append(linepack_savings(qin, qout, times, econ_data))
            LRL.append(linepack_loss(qout, mc, times, econ_data))
            MPR.append(market_power(qout, times, econ_data))
            CR.append(congestion_rent(qout, mc, times, econ_data))

            tks.append(times)
            prices.append([get_price(econ_data, tk, qoutk) for (tk, qoutk) in zip(times, qout)])
            costs.append([get_cost(econ_data, tk, qoutk) for (tk, qoutk) in zip(times, qout)])
            num_players = econ_data.number_players
            if num_players == 0:
                inv_num_players = 0
            else:
                inv_num_players = 1./num_players

            markups.append([-get_P_slope(econ_data, tk) * qoutk * inv_num_players for (tk, qoutk) in zip(times, qout)])

            # For there figure, we want to save the hdf5 associated with a data point

            if labels[0] == 'pmax':
                label2fnames[p_up_bar] = os.path.basename(fname)
            elif labels[0] == 'qlo':
                label2fnames[qlo] = os.path.basename(fname)

            fig_labels.append(pde_cournot.get_label(labels, cfg, args))

    if no_plot:
        return {'fnames': fnames, 'PS': PS, 'CS': CS, 'LRL': LRL, 'MPR': MPR, 'CR': CR, 'LCS': LCS,
                'pmax': pmax, 'qlo': qlos}

    d = datetime.now()
    date_suffix = d.isoformat(timespec='seconds').replace(':', '_')

    metadata = {'jsons': json.dumps(jsons_dict),
                'label2fnames': json.dumps(label2fnames),
                'cmd line': ' '.join(sys.argv)}

    if args.price_only:
        lLS = copy.deepcopy(LS)
        lcolor = copy.deepcopy(line_colors)

        fig_prices_lines, ax_prices_lines = plt.subplots()
        ax_prices_lines.plot(tks[0], costs[0], c=color_cost, label='Cost')

#        dat_prices = []
#        dat_prices.append({"vals": markups[idx], "color": color_markup, "label": "Markup"})
#        tk = tks[idx]
#        time_data = {'vals': np.asarray(tk), 'label': r'time'}

        for (idx, price) in enumerate(prices):
            ls = next(lLS)
            c = next(lcolor)

            line, = ax_prices_lines.plot(tks[idx], price, c=c, ls=ls, label=f'Price for $\\bar{{p}}$ = {pmax[idx]} bar', zorder=4)
            ax_prices_lines.plot(tks[idx], np.asarray(costs[idx])+np.asarray(markups[idx]), c=line.get_color(),
                    label=f'Static price for $\\bar{{p}}$ = {pmax[idx]} bar')

        ax_prices_lines.legend()
        fig_prices_lines.savefig(f"price_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)
        return

    xdict = {'pmax': pmax, 'qlo': qlos}
    xvals, xlabel = get_xaxis_data(labels, xdict)

    fig1, axs = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    ax_obj = axs[0, 0]
    ax_PS = axs[0, 1]
    ax_CS = axs[1, 0]
    ax_welfare = axs[1, 1]

    ax_obj.set_title("Objective value without penalization")
    ax_PS.set_title("Total producer surplus")
    ax_CS.set_title("Total consumer surplus")
    ax_welfare.set_title("Welfare")

    ax_obj.set_xlabel(xlabel)
    ax_PS.set_xlabel(xlabel)
    ax_CS.set_xlabel(xlabel)
    ax_welfare.set_xlabel(xlabel)

    fig_PS, axs = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    ax_LCS = axs[0, 0]
    ax_LRL = axs[0, 1]
    ax_MPR = axs[1, 0]
    ax_CR = axs[1, 1]

    ax_LCS.set_title("Linepack cost savings")
    ax_LRL.set_title("Linepack revenue loss")
    ax_MPR.set_title("Market power rent")
    ax_CR.set_title("Congestion rent")

    ax_obj.set_xlabel(xlabel)
    ax_LRL.set_xlabel(xlabel)
    ax_MPR.set_xlabel(xlabel)
    ax_CR.set_xlabel(xlabel)

    ax_obj.plot(xvals, obj, 'x-')
    ax_PS.plot(xvals, PS, 'x-')
    ax_CS.plot(xvals, CS, 'x-')
    ax_welfare.plot(xvals, welfare, 'x-')

    ax_LCS.plot(xvals, LCS, 'x-')
    ax_LRL.plot(xvals, LRL, 'x-')
    ax_MPR.plot(xvals, MPR, 'x-')
    ax_CR.plot(xvals, CR, 'x-')

    fig1.savefig(f"econ_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)
    fig_PS.savefig(f"PS_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    #fig_moor_econ, (ax_surplus, ax_prod_surplus) = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)

    dat1_pos = []
    dat1_neg = []
    dat2_pos = []

    dat1_pos.append({"vals": CR, "color": color_CR, "label": "Congestion rent"})
    dat1_pos.append({"vals": LCS, "color": color_LCS, "label": "Linepack savings"})
    dat1_pos.append({"vals": MPR, "color": color_MPR, "label": "Market power"})

    dat1_neg.append({"vals": LRL, "color": color_LRL, "label": "Linepack revenue loss"})

    dat2_pos.append({"vals": LCS, "color": color_LCS, "label": "Linepack savings"})
    dat2_pos.append({"vals": MPR, "color": color_MPR, "label": "Market power"})

    dat2_special = {"vals": np.asarray(LRL)+np.asarray(CR), "color": color_NR, "label": "Network rent"}

    dat_surplus = []
    dat_surplus.append({"vals": PS, "color": color_PS, "label": "Firms' profits"})
    dat_surplus.append({"vals": CS, "color": color_CS, "label": "Consumers surplus"})

    xdata = {'vals': np.asarray(xvals), 'label': xlabel}

    ############################################################################
    # Plot the various components of producer surplus
    ############################################################################

    if with_stack:
        fig_PS_stack1, ax1 = plt.subplots()
        plot_stack_with_total(ax1, xdata, dat1_pos, dat1_neg)
        fig_PS_stack1.savefig(f"PS_stack1_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

        fig_PS_stack2, ax2 = plt.subplots()
        plot_stack_with_total_and_special(ax2, xdata, dat2_pos, dat2_special)
        fig_PS_stack2.savefig(f"PS_stack2_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

        fig_welfare_stack, ax_welfare_stack = plt.subplots()
        plot_stack_with_total(ax_welfare_stack, xdata, dat_surplus, [], kw_total={'label': 'Welfare', 'color': color_W})
        fig_welfare_stack.savefig(f"welfare_stack_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    ############################################################################
    # Plot the same data with lines
    ############################################################################

    fig_PS1_lines, ax_PS1_lines = plt.subplots()
    fig_PS2_lines, ax_PS2_lines = plt.subplots()
    fig_PS_CS_welfare_lines, ax_PS_CS_welfare_lines = plt.subplots()
    fig_welfare_lines, ax_welfare_lines = plt.subplots()

    plot_lines_with_total(ax_PS1_lines, xdata, dat1_pos, dat1_neg, kw_total={'label': 'Firms\' profits', 'color': color_PS})

    plot_lines_with_total(ax_PS2_lines, xdata, dat2_pos, [dat2_special], kw_total={'label': 'Firms\' profits', 'color': color_PS})

    plot_lines_with_total(ax_PS_CS_welfare_lines, xdata, dat2_pos, [dat2_special], kw_total={'label': 'Firms\' profits', 'color': color_PS}, legend=False)
    ax_PS_CS_welfare_lines.plot(xvals, CS, c=color_CS, lw=lw_surplus, label='Consumers surplus')
    ax_PS_CS_welfare_lines.plot(xvals, np.asarray(PS)+np.asarray(CS), c=color_W, lw=lw_surplus, label='Welfare')

    plot_lines_with_total(ax_welfare_lines, xdata, dat_surplus, [], kw_total={'label': 'Welfare', 'color': color_W})

    fig_PS1_lines.savefig(f"PS1_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)
    fig_PS2_lines.savefig(f"PS2_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    ax_PS_CS_welfare_lines.legend()
    fig_PS_CS_welfare_lines.savefig(f"PS_CS_welfare_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    fig_welfare_lines.savefig(f"welfare_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    for (idx, xval) in enumerate(xvals):
        dat_prices = []
        dat_prices.append({"vals": costs[idx], "color": color_cost, "label": "Cost"})
        dat_prices.append({"vals": markups[idx], "color": color_markup, "label": "Markup"})
        tk = tks[idx]
        time_data = {'vals': np.asarray(tk), 'label': r'time'}

        fig_prices_stack, ax_prices_stack = plt.subplots()
        ax_prices_stack.plot(tk, prices[idx], c=color_price, label='Price', zorder=4)
        plot_stack_with_total(ax_prices_stack, time_data, dat_prices, [], kw_total={'label': 'Static price'})
        fig_prices_stack.savefig(f"price_stack_{date_suffix}_{'__'.join(labels)}-{xval}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

#        dat_prices.append({"vals": prices[idx], "color": 'r', "label": "Price"})
        fig_prices_lines, ax_prices_lines = plt.subplots()
        ax_prices_lines.plot(tk, prices[idx], c=color_price, label='Price', zorder=4)
        plot_lines_with_total(ax_prices_lines, time_data, dat_prices, [], kw_total={'label': 'Static price'})
        fig_prices_lines.savefig(f"price_lines_{date_suffix}_{'__'.join(labels)}-{xval}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)


def fig_dirs(args):

    results = {}
    for dir in args.dirs:
        hdf5_fnames = glob.glob(f"{dir}/*.hdf5")
        if not hdf5_fnames:
            print("No HDF5 file names given!")
            sys.exit(1)

        results[dir] = main(hdf5_fnames, args, no_plot=True)

    fig_PS1_lines, ax_PS1_lines = plt.subplots()
    fig_PS2_lines, ax_PS2_lines = plt.subplots()
    fig_PS_CS_welfare_lines, ax_PS_CS_welfare_lines = plt.subplots()
    fig_welfare_lines, ax_welfare_lines = plt.subplots()

    d = datetime.now()
    date_suffix = d.isoformat(timespec='seconds').replace(':', '_')

    # TODO: add these metadata
    # metadata = 'jsons': json.dumps(jsons_dict), 'label2fnames': json.dumps(label2fnames),
    metadata = {'cmd line': ' '.join(sys.argv)}

    dirs_len = len(args.dirs)
    allLSdirs = ("solid", "dashed", "dashdot", "dotted")
    LSdirs = dict(zip(args.dirs, islice(allLSdirs, dirs_len)))

    for dir in args.dirs:
        res = results[dir]
        m1 = re.match(r'.*N-([0-9]+)__.*', dir)
        m2 = re.match(r'.*-([0-9]+)-firms.*', dir)
        if 'welfare' in dir:
            dir_label = 'perf. comp.'
        if 'perf-comp' in dir:
            dir_label = 'perf. comp.'
        elif m1:
            dir_label = f'M = {m1.group(1)}'
        elif m2:
            dir_label = f'M = {m2.group(1)}'
        else:
            print(f"ERROR: unrecognized dir name {dir}. Must match 'perf-comp' or 'N-([0-9]+)__' or '-([0-9]+)-firms'")
            sys.exit(1)

        xdict = {'pmax': res['pmax'], 'qlo': res['qlo']}
        if not args.label:
            label = ['pmax']
        else:
            label = args.label

        xvals_unsorted, xlabel = get_xaxis_data(label, xdict)
        sidx = np.argsort(np.asarray(xvals_unsorted))
        xvals = np.asarray(xvals_unsorted)[sidx]

        dat1_pos = []
        dat1_neg = []
        dat2_pos = []

        CR = np.asarray(res['CR'])[sidx]
        LCS = np.asarray(res['LCS'])[sidx]
        MPR = np.asarray(res['MPR'])[sidx]
        LRL = np.asarray(res['LRL'])[sidx]
        PS = np.asarray(res['PS'])[sidx]
        CS = np.asarray(res['CS'])[sidx]

        dat1_pos.append({"vals": CR, "color": color_CR, "label": f"Network rent ({dir_label})"})
        dat1_pos.append({"vals": LCS, "color": color_LCS, "label": f"Linepack savings ({dir_label})"})
        dat1_pos.append({"vals": MPR, "color": color_MPR, "label": f"Market power ({dir_label})"})

        dat1_neg.append({"vals": LRL, "color": color_LRL, "label": f"Linepack revenue loss ({dir_label})"})

        dat2_pos.append({"vals": LCS, "color": color_LCS, "label": f"Linepack savings ({dir_label})"})
        dat2_pos.append({"vals": MPR, "color": color_MPR, "label": f"Market power ({dir_label})"})

        dat2_special = {"vals": np.asarray(LRL)+np.asarray(CR), "color": color_NR, "label": f"Network rent ({dir_label})"}

        dat_surplus = []
        dat_surplus.append({"vals": PS, "color": color_PS, "label": f"Firms' profits ({dir_label})"})
        dat_surplus.append({"vals": CS, "color": color_CS, "label": f"Consumers surplus ({dir_label})"})

        dat_CS = []
        dat_CS.append({"vals": CS, "color": color_CS, "label": f"Consumers surplus ({dir_label})"})

        xdata = {'vals': np.asarray(xvals), 'label': xlabel}

        kw_total = {'label': f"Firms' profits ({dir_label})", 'ls': LSdirs[dir], 'color': color_PS}
        kw_plots = {'ls': LSdirs[dir]}
        plot_lines_with_total(ax_PS1_lines, xdata, dat1_pos, dat1_neg, kw_plots=kw_plots, kw_total=kw_total)

        plot_lines_with_total(ax_PS2_lines, xdata, dat2_pos, [dat2_special], kw_plots=kw_plots, kw_total=kw_total)

        plot_lines_with_total(ax_PS_CS_welfare_lines, xdata, dat2_pos, [dat2_special], kw_plots=kw_plots, kw_total=kw_total, legend=False)
        ax_PS_CS_welfare_lines.plot(xvals, CS, c=color_CS, lw=lw_surplus, label="Consumers surplus")
        ax_PS_CS_welfare_lines.plot(xvals, np.asarray(PS)+np.asarray(CS), c=color_W, lw=lw_surplus, label='Welfare')

        ls = LSdirs[dir]
        ax_welfare_lines.plot(xvals, CS, c=color_CS, ls=ls, label=f"Consumers surplus ({dir_label})")
        ax_welfare_lines.plot(xvals, np.asarray(PS)+np.asarray(CS), c=color_W, ls=ls, label=f"Welfare ({dir_label})")

    fig_PS1_lines.savefig(f"PS1_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)
    fig_PS2_lines.savefig(f"PS2_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    ax_PS_CS_welfare_lines.legend()
    fig_PS_CS_welfare_lines.savefig(f"PS_CS_welfare_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)

    ax_welfare_lines.set_xlabel(xlabel)
    ax_welfare_lines.legend()
    fig_welfare_lines.savefig(f"welfare_lines_{date_suffix}.png", bbox_inches='tight', pad_inches=0, dpi=600, metadata=metadata)


if __name__ == '__main__':
    parser = pde_cournot.get_common_parser()
    parser.add_argument('--dirs', type=lambda s: [item for item in s.split(',')])
    parser.add_argument('--qout_price', type=lambda s: [item for item in s.split(',')])
    parser.add_argument('--price_only', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    hdf5_fnames = args.files

    if args.dirs:
        fig_dirs(args)
    else:
        if not hdf5_fnames:
            print("No HDF5 file names given!")
            sys.exit(1)

        main(hdf5_fnames, args)
