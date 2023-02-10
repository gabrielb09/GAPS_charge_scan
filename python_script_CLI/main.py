#!/usr/bin/env python

import argparse
import os as os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_config import *
from erf_function import *
from charge_scan_noinj import charge_scan_noinj
from charge_scan import charge_scan
from threshold_scan import threshold_scan
from compute_par_inj import get_parasitic_injection

parser = argparse.ArgumentParser(prog = 'chargescan', description = 'Analyzes Charge Scan Data from GAPS Tracker')
parser.add_argument('inFile', help = 'filepath to charge scan file')
parser.add_argument('-o', '--outPath', help = 'filepath to output directory')
parser.add_argument('-c', '--minChan', type = int, default = 0, help = 'minimum channel to analyze')
parser.add_argument('-C', '--maxChan', type = int, default = 31, help = 'maximum channel to analyze')
parser.add_argument('-e', '--exclude', type = int, nargs='*', help = 'channels to exclude')
parser.add_argument('-p', '--parasitic', action = 'store_true', help = 'compensate for parasitic capacitance')
parser.add_argument('-P', '--pedestal', help = 'filepath to pedestal data (necessary for parastitic capacitance correction)')
parser.add_argument('-T', '--transferFunc', help = 'filepath to transfer function data (necessary for parastitic capacitance correction)')
parser.add_argument('-t', '--peakingTime', type = int, help = 'peaking time')

args = parser.parse_args()

# LaTex interpreter
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})

# PLOT CONFIGURATION
# Label size
matplotlib.rcParams["axes.labelsize"] = 13
# Tick label size
matplotlib.rcParams["xtick.labelsize"] = 13
matplotlib.rcParams["ytick.labelsize"] = 13
# Figure size
matplotlib.rcParams["figure.figsize"] = 6.4 * 1.5, 4.8 * 1.5
# Legend font size
matplotlib.rcParams["legend.fontsize"] = 10

# Read data from file
data = pd.read_csv(
    args.inFile,
    comment="#",
    sep="\t",
    header=None,
)
data_bkp = data

# configure outpath
if args.outPath is None:
    args.outPath = args.inFile[:args.inFile.find('data')] + 'analysis'
    print(args.outPath)

# Configuration
conv_factor = 0.841
channels = range(args.minChan, args.maxChan + 1)
channels = np.setdiff1d(channels, args.exclude)

# Determine if charge scan or threshold scan
n_events = data.iloc[0][2]
threshold_col = data.iloc[:, 0].to_numpy()
thr_unique = np.unique(threshold_col)
charge_scan_flag = True

if len(thr_unique) > 1:
    charge_scan_flag = False

if charge_scan_flag:
    # Charge scan without removal of parasitic injection
    # Always done
    (xmin, xmax) = charge_scan(data, channels, conv_factor, args.outPath)

    if args.parasitic:
        if args.pedestal is None or args.transferFunc is None or args.peakingTime is None:
            raise(FileExistsError('Need pedestal, and transfer function data for parasitic injection compensation'))
        allch_par_inj_estimate = []
        for ch in channels:
            ch_par_inj_estimate = get_parasitic_injection(args.pedestal, args.transferFunc, ch, args.peakingTime)
            allch_par_inj_estimate.append(ch_par_inj_estimate)

        # Charge scan with subtracted parasitic injection
        charge_scan_noinj(
            data,
            channels,
            conv_factor,
            args.outPath,
            xmin,
            xmax,
        )

else: # Threshold scan
    threshold_scan(data_bkp, channels, n_events, output_folder_filepath)
