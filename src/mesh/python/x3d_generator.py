#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/x3d_generator.py
# date  Monday, Jul 19, 2021, 12:14 pm
# brief This script generates X3D mesh files from a mesh object (assumed to have certain data).
# note  Copyright (C) 2021, Triad National Security, LLC.,  All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import mesh_types
import numpy as np
import argparse

# -- mesh class dictionary
mesh_type_dict = {'orth_2d_mesh': mesh_types.orth_2d_mesh}

# ------------------------------------------------------------------------------------------------ #
# -- create argument parser

parser = argparse.ArgumentParser(description='Generate X3D mesh file.')
parser.add_argument('-mt', '--mesh_type', type=str, default='orth_2d_mesh',
                    help='Select mesh type.')
parser.add_argument('-nd', '--num_per_dim', type=int, nargs='+', default=[1, 1],
                    help='Number of cells per dimension.')
parser.add_argument('-bd', '--bnd_per_dim', type=float, nargs='+', default=[0.0, 1.0, 0.0, 1.0],
                    help='Length per dimension.')

# -- parse arguments from command line
args = parser.parse_args()

# ------------------------------------------------------------------------------------------------ #
# -- create the relevant mesh object

# -- number of spatial dimensions should be length of num_per_dim
ndim = len(args.num_per_dim)

# -- sanity check input
assert (ndim > 0), 'len(args.num_per_dim) <= 0'
assert (ndim < 4), 'len(args.num_per_dim) >= 4'
assert (len(args.bnd_per_dim) == 2 * ndim), 'len(args.bnd_per_dim) != 2 * ndim'

# -- de-serialize spatial bound list
bnd_per_dim = [[args.bnd_per_dim[2 * i], args.bnd_per_dim[2 * i + 1]] for i in range(ndim)]

# -- instantiate the class for the mesh type selected
mesh = mesh_type_dict[args.mesh_type](bnd_per_dim, args.num_per_dim)

# ------------------------------------------------------------------------------------------------ #
# -- write out the main mesh file in x3d format

# -- open writable file
fo = open('x3d.mesh.in', 'w')

# -- write header block
# -- for x3d:
# --  nodes_per_face = max number of nodes per face
# --  faces_per_cell = max number of faces per cell
fo.write('ascii\n')
fo.write('header\n')
fo.write('   process                         1\n')
fo.write('   numdim                          {0}\n'.format(mesh.ndim))
fo.write('   materials                       1\n')
fo.write('   nodes                           {0}\n'.format(mesh.num_nodes))
fo.write('   faces                           {0}\n'.format(mesh.num_faces))
fo.write('   elements                        {0}\n'.format(mesh.num_cells))
fo.write('   ghost_nodes                     0\n')
fo.write('   slaved_nodes                    0\n')
fo.write('   nodes_per_slave                 2\n')
fo.write('   nodes_per_face                  {0}\n'.format(np.max(mesh.num_nodes_per_face)))
fo.write('   faces_per_cell                  {0}\n'.format(np.max(mesh.num_faces_per_cell)))
fo.write('   node_data_fields                0\n')
fo.write('   cell_data_fields                2\n')
fo.write('end_header\n')

# -- space
fo.write('\n')

# -- write trivial material block
fo.write('matnames\n')
fo.write('            1   0\n')
fo.write('end_matnames\n')

# -- space
fo.write('\n')

# -- write trivial material eos block
fo.write('mateos\n')
fo.write('            1   -1\n')
fo.write('end_mateos\n')

# -- space
fo.write('\n')

# -- write trival material opc block
fo.write('matopc\n')
fo.write('            1   -1\n')
fo.write('end_matopc\n')

# -- space
fo.write('\n')

# -- write node coordinate block
fo.write('nodes\n')
for node in range(mesh.num_nodes):
    crds = [0.0, 0.0, 0.0]
    for idim in range(mesh.ndim):
        crds[idim] = mesh.coordinates_per_node[node, idim]
    fo.write('{0:10d}  {1:.15e}  {2:.15e}  {3:.15e}\n'.format(node + 1, crds[0], crds[1], crds[2]))
fo.write('end_nodes\n')

# -- space
fo.write('\n')

# -- write face block (excludes extra columns used for parallel meshes, for now)
fo.write('faces\n')
for cell in range(mesh.num_cells):
    for j in range(mesh.num_faces_per_cell[cell]):
        # -- get the full face index for the cell and local face index
        face = mesh.faces_per_cell[cell, j]
        # -- get the number of nodes for this face
        nnpf = mesh.num_nodes_per_face[face]
        # -- build face string
        face_str = '{0:10d}{1:10d}'.format(face + 1, nnpf)
        # -- append node indices to face string
        for node in mesh.nodes_per_face[face, :]:
            face_str += '{0:10d}'.format(node + 1)
        # -- write face line
        fo.write(face_str + '\n')
fo.write('end_faces\n')

# -- space
fo.write('\n')

# -- write cell block
fo.write('cells\n')
for cell in range(mesh.num_cells):
    # -- get the number of faces for this cell
    nfpc = mesh.num_faces_per_cell[cell]
    # -- build cell string
    cell_str = '{0:10d}{1:10d}'.format(cell + 1, nfpc)
    # -- append face indices to cell string
    for face in mesh.faces_per_cell[cell, :]:
        cell_str += '{0:10d}'.format(face + 1)
    # -- write cell line
    fo.write(cell_str + '\n')
fo.write('end_cells\n')

# -- space
fo.write('\n')

# -- unused AMR-type/parallel data
fo.write('slaved_nodes         0\n')
fo.write('end_slaved_nodes\n')
fo.write('ghost_nodes          0\n')
fo.write('end_ghost_nodes\n')

# -- space
fo.write('\n')

# -- unused data per cell
fo.write('cell_data\n')
fo.write('matid\n')
fo.write('         0\n')
fo.write('end_matid\n')
fo.write('partelm\n')
fo.write('         1\n')
fo.write('end_partelm\n')
fo.write('end_cell_data\n')

# -- close main mesh file
fo.close()

# ------------------------------------------------------------------------------------------------ #
# -- write out the mesh boundary node lists to file (for boundary conditions)

# -- assume two boundaries per dimension
# -- x3d only needs a unique node list
for i in range(2 * mesh.ndim):
    np.savetxt('x3d.mesh.bdy' + str(i + 1) + '.in', np.unique(mesh.nodes_per_side[i] + 1), fmt='%d')

# ------------------------------------------------------------------------------------------------ #
# end of x3d_generator.py
# ------------------------------------------------------------------------------------------------ #
