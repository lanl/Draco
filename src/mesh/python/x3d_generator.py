#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/x3d_generator.py
# date  Monday, Jul 19, 2021, 12:14 pm
# brief This script generates X3D mesh files from a mesh object (assumed to have certain data).
# note  Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import mesh_types
import numpy as np
import argparse

# -- mesh class dictionary
mesh_type_dict = {'orth_1d_mesh': mesh_types.orth_1d_mesh, 'orth_2d_mesh': mesh_types.orth_2d_mesh,
                  'orth_3d_mesh': mesh_types.orth_3d_mesh, 'fcc_3d_mesh': mesh_types.fcc_3d_mesh,
                  'vor_2d_mesh': mesh_types.vor_2d_mesh, 'rnd_1d_mesh': mesh_types.rnd_1d_mesh,
                  'rnd_2d_mesh': mesh_types.rnd_2d_mesh, 'rnd_3d_mesh': mesh_types.rnd_3d_mesh}

# ------------------------------------------------------------------------------------------------ #
# -- create argument parser

parser = argparse.ArgumentParser(description='Generate X3D mesh file.')
parser.add_argument('-mt', '--mesh_type', type=str, default='orth_2d_mesh',
                    help='Select mesh type.')
parser.add_argument('-nd', '--num_per_dim', type=int, nargs='+', default=[1, 1],
                    help='Number of cells per dimension.')
parser.add_argument('-bd', '--bnd_per_dim', type=float, nargs='+', default=[0.0, 1.0, 0.0, 1.0],
                    help='Bounds per dimension.')
parser.add_argument('--name', type=str, default='mesh',
                    help='Select file name (will be prefixed with x3d. and sufficed with .in).')
parser.add_argument('--num_cells', type=int, default=100,
                    help='Number of cells')
# -- random mesh options
parser.add_argument('--eps', type=float, default=1.0,
                    help='0=no randomization, 1=max possible randomization.')
parser.add_argument('--rnd_seed', type=int, default=1,
                    help='Random number seed (for random mesh).')
# -- optional region/matid input
parser.add_argument('--reg_ids', type=int, nargs='+', default=[], help='List of region ids.')
parser.add_argument('--reg_bnd_per_dim', type=float, nargs='+', default=[],
                    help='Bounds per dimension of region.')

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

# -- sanity check optional region input
assert (len(args.reg_bnd_per_dim) == 2 * ndim * len(args.reg_ids)), \
    'len(args.reg_bnd_per_dim) != 2 * ndim * len(args.reg_ids)'

# -- de-serialize spatial bound list
bnd_per_dim = [[args.bnd_per_dim[2 * i], args.bnd_per_dim[2 * i + 1]] for i in range(ndim)]

# -- instantiate the class for the mesh type selected
if args.mesh_type in ['orth_1d_mesh', 'orth_2d_mesh', 'orth_3d_mesh', 'fcc_3d_mesh']:
    mesh = mesh_type_dict[args.mesh_type](bnd_per_dim, args.num_per_dim)
elif args.mesh_type in ['vor_2d_mesh']:
    mesh = mesh_type_dict[args.mesh_type](bnd_per_dim, args.num_cells)
elif args.mesh_type in ['rnd_1d_mesh', 'rnd_2d_mesh', 'rnd_3d_mesh']:
    mesh = mesh_type_dict[args.mesh_type](bnd_per_dim, args.num_per_dim, args.eps, args.rnd_seed)

# ------------------------------------------------------------------------------------------------ #
# -- construct matid and region dictionary if region data arguments are available

# -- set matid data (using region ids and bounds if available)
matids = np.zeros(mesh.num_cells, dtype=int)
reg_dict = {}
reg_ids = []
if (len(args.reg_ids) > 0):
    # -- prepend reg_ids and reg_bnd_per_dim with universal region
    reg_ids = args.reg_ids
    reg_bnd_per_dim = args.reg_bnd_per_dim
    reg_ids = [0] + reg_ids
    reg_bnd_per_dim = args.bnd_per_dim + reg_bnd_per_dim
    # -- initialize dictionary of cell lists
    reg_dict = {reg_id: [] for reg_id in reg_ids}
    nmats = len(reg_ids)
    for cell in range(mesh.num_cells):
        # TODO: make coordinate averaging a method for each mesh_type (permit distinct approaches)
        # -- initialize container for cell average coordinate
        coord_av = np.array([0.0, 0.0, 0.0])
        # -- average coordinates over faces first
        for face in mesh.faces_per_cell[cell]:
            # -- initialize container for face average coordinate
            coord_face_av = np.array([0.0, 0.0, 0.0])
            # -- get the node coordinates for this face to average
            for node in mesh.nodes_per_face[face]:
                for idim in range(mesh.ndim):
                    coord_face_av[idim] += mesh.coordinates_per_node[node, idim]
            # -- finish face average by dividing by number of nodes per face
            coord_face_av /= mesh.num_nodes_per_face[face]
            # -- add face average to cell average
            coord_av += coord_face_av
        # -- finish cell average by dividing by number of faces per cell
        coord_av /= mesh.num_faces_per_cell[cell]
        # -- check which region to which this cell belongs
        for imat in range(nmats):
            mat_has_cell = True
            # -- check coord_av against region bounds along each dimension
            m_off = 2 * ndim * imat
            for idim in range(mesh.ndim):
                mat_has_cell = mat_has_cell and \
                    coord_av[idim] >= reg_bnd_per_dim[2 * idim + m_off] and \
                    coord_av[idim] < reg_bnd_per_dim[2 * idim + 1 + m_off]
            # -- set mat/reg data if the cell is in the region
            if (mat_has_cell):
                matids[cell] = reg_ids[imat]
                reg_dict[reg_ids[imat]].append(cell + 1)

    # -- remove region cells from universal region and lower-index regions
    # -- this implies that this region takes precedence over regions earlier in the iteration
    for j in range(nmats - 1):
        for reg_id in reg_ids[j + 1:]:
            reg_dict[j] = [cell for cell in reg_dict[j] if cell not in reg_dict[reg_id]]

    # -- convert dictionary of lists to dictionary of numpy arrays
    for reg_id in reg_ids:
        reg_dict[reg_id] = np.array(reg_dict[reg_id], dtype=int)

else:
    reg_ids = [0]
    reg_dict[0] = np.arange(1, mesh.num_cells + 1)

# ------------------------------------------------------------------------------------------------ #
# -- write out the main mesh file in x3d format

# -- open writable file
fname = 'x3d.' + args.name
fo = open(fname + '.in', 'w')

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
        face = mesh.faces_per_cell[cell][j]
        # -- get the number of nodes for this face
        nnpf = mesh.num_nodes_per_face[face]
        # -- build face string
        face_str = '{0:10d}{1:10d}'.format(face + 1, nnpf)
        # -- append node indices to face string
        for node in mesh.nodes_per_face[face]:
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
    for face in mesh.faces_per_cell[cell]:
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

# -- write matids per cell
fo.write('cell_data\n')
fo.write('matid\n')
for cell in range(mesh.num_cells):
    fo.write('{0:10d}{1:10d}\n'.format(cell + 1, matids[cell]))
fo.write('end_matid\n')

# -- unused data (possibly supposed to be per cell)
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
    np.savetxt(fname + '.bdy' + str(i + 1) + '.in',
               np.unique(np.array(mesh.nodes_per_side[i]) + 1), fmt='%d')

# ------------------------------------------------------------------------------------------------ #
# -- write out the mesh region files based on matid
# -- todo: map (dictionary) matid to cell lists (this is only for X3D, so it's not a mesh task)
j = 0
for i in range(len(reg_ids)):
    # -- only save if region has cells
    if (np.size(reg_dict[reg_ids[i]]) > 0):
        np.savetxt(fname + '.Reg' + str(j + 1) + '.in', reg_dict[reg_ids[i]], fmt='%d')
        # -- incrment region counter
        j += 1

# ------------------------------------------------------------------------------------------------ #
# end of x3d_generator.py
# ------------------------------------------------------------------------------------------------ #
