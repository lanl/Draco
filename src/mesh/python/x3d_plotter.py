#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/x3d_plotter.py
# date  Monday, Aug 9, 2021
# brief This script plots X3D mesh files.
# note  Copyright (C) 2021, Triad National Security, LLC.,  All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import matplotlib.pyplot as plt
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

# Data to read in
numdim = None
numnodes = None
numfaces = None
numcells = None
nodes = []
faces = []
cells = []

blocks = ['header', 'nodes', 'faces', 'cells']
current_block = None
for line in lines:
    words = line.split()
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
    if current_block == 'header':
        if words[0] == 'numdim':
            numdim = int(words[1])
        elif words[0] == 'nodes':
            numnodes = int(words[1])
        elif words[0] == 'faces':
            numfaces = int(words[1])
        elif words[0] == 'elements':
            numcells = int(words[1])
    elif current_block == 'nodes':
        if len(words) == 4:
            nodes.append([float(words[1]), float(words[2]), float(words[3])])
    elif current_block == 'faces':
        if len(words) >= 3:
            face = []
            for nnodes in range(int(words[1])):
                face.append(int(words[nnodes + 2]) - 1)
            faces.append(face)
    elif current_block == 'cells':
        if len(words) >= 3:
            cell = []
            for nface in range(int(words[1])):
                cell.append(int(words[nface + 2]) - 1)
            cells.append(cell)

assert (numdim is not None), "numdim not found!"
assert (numnodes is not None), "numnodes not found!"
assert (numfaces is not None), "numfaces not found!"
assert (numcells is not None), "numcells not found!"
assert (len(nodes) == numnodes), "numnodes does not match number of nodes!"
assert (len(faces) == numfaces), "numfaces does not match number of faces!"
assert (len(cells) == numcells), "numcells does not match number of faces!"

# ------------------------------------------------------------------------------------------------ #
# -- Plot mesh

plt.figure()
ax = plt.gca()

for cell in cells:
    for face in cell:
        pt1 = nodes[faces[face][0]]
        pt2 = nodes[faces[face][1]]
        ax.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], color='k')
plt.show()
