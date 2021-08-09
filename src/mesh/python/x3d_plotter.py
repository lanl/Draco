#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/x3d_plotter.py
# date  Monday, Aug 9, 2021
# brief This script plots X3D mesh files.
# note  Copyright (C) 2021, Triad National Security, LLC.,  All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import mesh_types
import numpy as np
import argparse
import os

# ------------------------------------------------------------------------------------------------ #
# -- create argument parser

parser = argparse.ArgumentParser(description='Plot X3D mesh file.')
parser.add_argument('-fn', '--file_name', type=str, default=None, required=True,
                    help='Provide mesh file to plot.')

# -- parse arguments from command line
args = parser.parse_args()

# ------------------------------------------------------------------------------------------------ #
# -- Read and parse x3d file

assert (os.path.exists(args.file_name)), f"Mesh file \"{args.file_name}\" does not exist!"
with open(args.file_name) as f:
    lines = [line.strip() for line in f]

print(lines)

blocks = ['header', 'nodes', 'faces', 'cells']
current_block = None
for line in lines:
    print(line)
    # If no current block, check if starting new block
    if current_block is None:
        for block in blocks:
            if block == line:
                current_block = block
                break
    # If current block, check if ending current block
    else:
        if line == "end_" + current_block:
            current_block = None

    # Process data if currently on a block
    if current_block is not None:
        print(current_block)
