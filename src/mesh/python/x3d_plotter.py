#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/x3d_plotter.py
# date  Monday, Aug 9, 2021
# brief This script plots X3D mesh files.
# note  Copyright (C) 2021, Triad National Security, LLC.,  All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import matplotlib.pyplot as plt
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
face_indices = []
faces = []
cells = []
boundaries = []

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
                # Convert from file node ID to code node index
                face.append(int(words[nnodes + 2]) - 1)
            face_index = int(words[0])
            # if face_index not in face_indices:
            faces.append(face)
            face_indices.append(int(words[0]))
    elif current_block == 'cells':
        if len(words) >= 3:
            cell = []
            for nface in range(int(words[1])):
                # Convert from file face ID to code face index
                cell.append(int(words[nface + 2]) - 1)
            cells.append(cell)

# Sort faces in case they are out of order
faces = [x for _, x in sorted(zip(face_indices, faces))]

# Read boundaries
boundary_files = []
boundary_nodes = []
boundary_faces = []
if numdim == 2:
    for n in range(4):
        assert (args.file_name[-3:] == '.in'), "Filename does not end in \".in\""
        boundary_files.append(args.file_name[:-3] + f".bdy{n+1}.in")

for boundary_file in boundary_files:
    with open(boundary_file) as f:
        lines = [line.strip() for line in f]
    # -- read in boundary nodes
    boundary = []
    for line in lines:
        boundary.append(int(line) - 1)
    boundary_nodes.append(boundary)

    # -- calculate boundary faces
    boundary_face_tmp = []
    for face_idx, face in enumerate(faces):
        node0 = face[0]
        node1 = face[1]
        if node0 in boundary and node1 in boundary:
            boundary_face_tmp.append(face_idx)
    boundary_faces.append(boundary_face_tmp)


# -- sanity checks
assert (numdim is not None), "numdim not found!"
assert (numnodes is not None), "numnodes not found!"
assert (numfaces is not None), "numfaces not found!"
assert (numcells is not None), "numcells not found!"
assert (len(nodes) == numnodes), "numnodes does not match number of nodes!"
assert (len(faces) == numfaces), "numfaces does not match number of faces!"
assert (len(cells) == numcells), "numcells does not match number of faces!"

# ------------------------------------------------------------------------------------------------ #
# -- Plot mesh

if numdim == 1:
    assert (False), "1D plotting not supported!"
elif numdim == 2:

    plt.figure()
    ax = plt.gca()

    # -- plot faces
    plotted_faces = []
    for cell in cells:
        for face in cell:
            # Don't re-plot the same face
            if (([faces[face][0], faces[face][1]] not in plotted_faces) and
                ([faces[face][1], faces[face][0]] not in plotted_faces)):
                pt1 = nodes[faces[face][0]]
                pt2 = nodes[faces[face][1]]
                plotted_faces.append([faces[face][0], faces[face][1]])
                ax.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], color='k')

    # -- plot boundary faces
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    for n, bound in enumerate(boundary_faces):
        for face in bound:
            pt1 = nodes[faces[face][0]]
            pt2 = nodes[faces[face][1]]
            ax.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], color=colors[n], linewidth=4)

    # -- plot nodes
    for node in nodes:
        ax.plot([node[0]], [node[1]], marker='.', color='b')
    plt.show()

elif numdim == 3:
    assert (False), "3D plotting not supported!"
