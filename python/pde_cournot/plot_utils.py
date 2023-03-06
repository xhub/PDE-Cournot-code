# SPDX-License-Identifier: MIT
import argparse
import matplotlib


def get_label_pmax(cfg: dict[str, dict], args) -> str:
    pipe = cfg['pipe_json']
    pmax_str = args.pressure_format.format(pipe.pressure_max/1e5)
    return r'$\bar{{p}}$ = {:} bar'.format(pmax_str)


def get_label_pinit(cfg: dict[str, dict], args) -> str:
    approx = cfg['approx_json']
    pinit_str = args.pressure_format.format(approx.operational_pressure/1e5)
    return r'$p_0(0)$ = {:} bar'.format(pinit_str)


def get_label_cost_diff(cfg: dict[str, dict], args) -> str:
    econ = cfg['econ_json']
    cost1 = econ.cost_functions[0]["intercept"]
    cost2 = econ.cost_functions[1]["intercept"]
    δcost = (cost2-cost1)/cost2

    cost_str = args.deltacost_format.format(δcost)

    return r'$\delta cost$ = {:} %'.format(cost_str)


def get_label_qlo(cfg: dict[str, dict], args) -> str:
    pipe = cfg['pipe_json']
    qlo = pipe.flow_min
    return r'$\underaccent{{\bar}}{{q}}$ = {:}'.format(qlo)


label2fn = {"pmax": get_label_pmax,
            "pinit": get_label_pinit,
            "cost_diff": get_label_cost_diff,
            "qlo": get_label_qlo,
            }


def get_label(labels: list[str], cfg: dict[str, dict], args) -> str:
    return ', '.join((label2fn[label](cfg, args) for label in labels))


def get_common_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--label', action='append')
    parser.add_argument('files', nargs='*')
    parser.add_argument('--colors', type=lambda s: [item for item in s.split(',')])
    parser.add_argument('--pressure-format', default='{:.1f}')
    parser.add_argument('--deltacost-format', default='{:.1f}')

    return parser


def set_latex_preamble():
    matplotlib.rcParams['text.usetex'] = True
    matplotlib.rcParams['text.latex.preamble'] = r'''\usepackage{amsmath}
        \usepackage{amssymb}
        \usepackage{accents}
'''
