# -*- coding: utf-8 -*-
"""
function is used to translate rinex into snr file format
These files are then used by gnss_lomb.py
author: kristine larson
date: 20 march 2019
"""
import sys
import os
import numpy as np
import gps
import argparse
import datetime

# set an environment variable for where you are keeping your LSP
# instructions and input files 
# DEFINE REFL_CODE on your system
#xdir = os.environ['REFL_CODE']
# this should be done in your .bashrc file

#WARNING - THIS CODE ASSUMES IF YOU USE HATANAKA RINEX
# FILES, YOU NEED THE APPROPRIATE EXECUTABLE
# SEE rinex_unavco.py for details
 
#
# user inputs the observation file information
parser = argparse.ArgumentParser()
parser.add_argument("station", help="station name", type=str)
parser.add_argument("year", help="year", type=int)
parser.add_argument("doy1", help="doy1", type=int)
parser.add_argument("doy2", help="doy2", type=int)
parser.add_argument("snrEnd", help="snrEnd", type=str)
parser.add_argument("orbType", help="orbType", type=str)
args = parser.parse_args()
#
# rename the user inputs as variables
#
station = args.station
year = args.year
doy1= args.doy1
doy2= args.doy2
snrt = args.snrEnd
orbtype = args.orbType
#

doy_list = list(range(doy1, doy2+1))

# for each day in the doy list
for doy in doy_list:
    gps.quick_rinex_snr(year, doy, station, snrt, orbtype)
