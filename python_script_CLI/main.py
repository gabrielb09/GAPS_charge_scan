#!/usr/bin/env python

import argparse
import os as os
import numpy as np
import pandas as pd
from pathlib import Path as path

from plot_config import *
from erf_function import *
from charge_scan_noinj import charge_scan_noinj
from charge_scan import charge_scan
from threshold_scan import threshold_scan
from compute_par_inj import get_parasitic_injection

parser = argparse.ArgumentParser(prog = 'chargescan', description = 'Analyzes Charge Scan Data from GAPS Tracker')
parser.add_argument('inFile', help = 'filepath to charge scan file')
parser.add_argument('outPath', help = 'filepath to output directory')
parser.add_argument('-c', '--minChan', type = int, default = 0, help = 'minimum channel to analyze')
parser.add_argument('-C', '--maxChan', type = int, default = 31, help = 'maximum channel to analyze')
parser.add_argument('-e', '--exclude', type = int, nargs='*', help = 'channels to exclude')
parser.add_argument('-p', '--parasitic', action = 'store_true', help = 'compensate for parasitic capacitance')
parser.add_argument('-P', '--pedestal', help = 'filepath to pedestal data (necessary for parastitic capacitance correction)')
parser.add_argument('-T', '--transferFunc', help = 'filepath to transfer function data (necessary for parastitic capacitance correction)')
parser.add_argument('-t', '--peakingTime', type = int, help = 'peaking time')

<<<<<<< HEAD
while True:
    print("\n*** GAPS CHARGE SCAN TOOL v1.2 ***\n")

    # Request user input
    filename_chargescan_path_flag = False
    while not filename_chargescan_path_flag:
        filename_chargescan = input("    Charge or threshold scan filepath: ")
        if filename_chargescan[0] == '"':
            filename_chargescan = filename_chargescan.replace('"', "")
        # Check if file exists
        filename_chargescan_path_flag = path(filename_chargescan).is_file()
        if not filename_chargescan_path_flag:
            print("\nInvalid filepath!\n")

    chmin_check_flag = False
    while not chmin_check_flag:
        ch_min = int(input("                        First channel: "))
        chmin_check_flag = ch_min >= 0 and ch_min <= 31
        if not chmin_check_flag:
            print("\nInvalid channel number!\n")

    chmax_check_flag = False
    while not chmax_check_flag:
        ch_max = int(input("                         Last channel: "))
        chmax_check_flag = ch_max >= 0 and ch_max <= 31
        if not chmax_check_flag:
            print("\nInvalid channel number!\n")

    excl_channels_check_flag = True
    while excl_channels_check_flag:
        excl_channels_check_flag = False
        excl_channels = input("  Excluded channels (comma separated): ")

        if excl_channels != "":
            excl_channels = excl_channels.split(",")
            excl_channels = [int(i) for i in excl_channels]

            for ch in excl_channels:
                temp_flag = ch < ch_min or ch > ch_max
                excl_channels_check_flag = excl_channels_check_flag or temp_flag

            if excl_channels_check_flag:
                print("\nInvalid channel numbers in list!\n")

        else:
            excl_channels = []
=======
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
>>>>>>> fcef162 (changed user input method to argparse)

# Read data from file
data = pd.read_csv(
    args.inFile,
    comment="#",
    sep="\t",
    header=None,
)
data_bkp = data

<<<<<<< HEAD
    # Ask threshold for channel deactivation
    deactivate_thr = input("      Deactivate channels below [keV]: ")
    if deactivate_thr != "":
        deactivate_thr = int(deactivate_thr)
    else:
        deactivate_thr = -500

    # Ask ENC for channel deactivation
    deactivate_enc = input(" Deactivate channels with ENC > [keV]: ")
    if deactivate_enc != "":
        deactivate_enc = int(deactivate_enc)
    else:
        deactivate_enc = 500
=======
# Configuration
conv_factor = 0.841
channels = range(args.minChan, args.maxChan + 1)
channels = np.setdiff1d(channels, args.exclude)

# Determine if charge scan or threshold scan
n_events = data.iloc[0][2]
threshold_col = data.iloc[:, 0].to_numpy()
thr_unique = np.unique(threshold_col)
charge_scan_flag = True
>>>>>>> fcef162 (changed user input method to argparse)

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

<<<<<<< HEAD
    if charge_scan_flag:
        # Ask user to compensate parasitic injection or not
        comp_parinj_check = False
        while not comp_parinj_check:
            comp_parinj_flag = input("Compensate parasitic injection? (y/n): ")
            comp_parinj_check = comp_parinj_flag == "y" or comp_parinj_flag == "n"
            if not comp_parinj_check:
                print('\nInvalid choice! Only "y" or "n" allowed!\n')

        # Charge scan with removal of parasitic injection
        if comp_parinj_flag == "y":
            # Get additional info when charge scan is selected and user wants to compensate parasitic injection
            pedestal_check_flag = False
            while not pedestal_check_flag:
                filename_pedestal = input("         Pedestal from automated test: ")
                if filename_pedestal[0] == '"':
                    filename_pedestal = filename_pedestal.replace('"', "")
                pedestal_check_flag = path(filename_pedestal).is_file()
                if not pedestal_check_flag:
                    print("\nInvalid filepath!\n")

            fdt_check_flag = False
            while not fdt_check_flag:
                filename_fdt = input("Transfer function from automated test: ")
                if filename_fdt[0] == '"':
                    filename_fdt = filename_fdt.replace('"', "")
                fdt_check_flag = path(filename_fdt).is_file()
                if not fdt_check_flag:
                    print("\nInvalid filepath!\n")

            pt_check_flag = False
            while not pt_check_flag:
                peaking_time = int(input("                Peaking time (0 to 7): "))
                pt_check_flag = peaking_time >= 0 and peaking_time <= 7
                if not pt_check_flag:
                    print("\nInvalid peaking time!\n")

        # Charge scan without removal of parasitic injection
        # Always done
        (xmin, xmax) = charge_scan(
            data,
            channels,
            conv_factor,
            output_folder_filepath,
            deactivate_thr,
            deactivate_enc,
        )

        # Charge scan with removal of parasitic injection
        if comp_parinj_flag == "y":
            # Write estimated parameters to file
            summary_data_filepath = os.path.join(
                output_folder_filepath,
                "summary_inj_ch" + str(ch_min) + "-" + str(ch_max) + ".dat",
            )
            with open(summary_data_filepath, "w") as fp:
                fp.write("ch\ttf_gain\ttf_pedestal\tauto_pedestal\tpar_inj\n")
            fp.close()

            # Get parasitic injection estimate to get proper charge scan
            allch_par_inj_estimate = []
            for ch in channels:
                ch_par_inj_estimate = get_parasitic_injection(
                    filename_pedestal,
                    filename_fdt,
                    ch,
                    peaking_time,
                    summary_data_filepath,
                )
                allch_par_inj_estimate.append(ch_par_inj_estimate)

            # Charge scan with subtracted parasitic injection
            charge_scan_noinj(
                data,
                channels,
                conv_factor,
                output_folder_filepath,
                xmin,
                xmax,
                allch_par_inj_estimate,
                deactivate_thr,
                deactivate_enc,
            )

    else:
        # Threshold scan
        threshold_scan(data_bkp, channels, n_events, output_folder_filepath)
=======
else: # Threshold scan
    threshold_scan(data_bkp, channels, n_events, output_folder_filepath)
>>>>>>> fcef162 (changed user input method to argparse)
