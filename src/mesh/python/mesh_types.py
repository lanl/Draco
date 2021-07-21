#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file  src/mesh/python/mesh_types.py
# date  Monday, Jul 19, 2021, 12:14 pm
# brief This script provides mesh classes that calculate and contain unstructred mesh data.
# note  Copyright (C) 2021, Triad National Security, LLC.,  All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import numpy as np


# ------------------------------------------------------------------------------------------------ #
# base class (to show required member data)
class base_mesh:
    '''
    Base class with data required by mesh file generators.
    This contains no methods for creating the data members.
    '''
    def __init__(self):
        # -- required data
        self.ndim = 0  # number of dimensions
        self.num_nodes = 0  # total number of nodes
        self.coordinates_per_node = np.array([])  # coordinate array indexed by node
        self.num_cells = 0  # total number of cells
        self.num_faces = 0  # total number of oriented faces
        self.num_faces_per_cell = np.array([], dtype=int)  # number of faces per cell
        self.num_nodes_per_face = np.array([], dtype=int)  # number of nodes per face
        self.faces_per_cell = np.array([[]], dtype=int)  # face indexes per cell
        self.nodes_per_face = np.array([[]], dtype=int)  # node indexes per face
        self.nodes_per_side = [np.array([[]], dtype=int)]  # list of arrays of node per bdy face


# ------------------------------------------------------------------------------------------------ #
# orthogonal 2D mesh type
class orth_2d_mesh(base_mesh):
    '''
    Class for orthogonally structured 2D mesh data.
    This class generates an orthogonally structured mesh in an unstructured
    format suitable for creating unstructured-mesh input files.
    '''
    def __init__(self, bounds_per_dim, num_cells_per_dim):

        # -- short-cuts
        nx = num_cells_per_dim[0]
        ny = num_cells_per_dim[1]

        # -- number of dimensions
        ndim = len(num_cells_per_dim)
        self.ndim = ndim
        assert (ndim == 2), 'ndim != 2, exiting...'
        assert (len(bounds_per_dim) == ndim), 'len(bounds_per_dim) != ndim, exiting...'

        # -- create grid arrays along each dimension
        grid_per_dim = [np.linspace(bounds_per_dim[i][0], bounds_per_dim[i][1],
                                    num_cells_per_dim[i] + 1) for i in range(ndim)]

        # -- create node indices
        num_nodes = (nx + 1) * (ny + 1)
        self.num_nodes = num_nodes
        self.coordinates_per_node = np.zeros((num_nodes, ndim))
        for j in range(ny + 1):
            j_offset = (nx + 1) * j
            for i in range(nx + 1):
                node = i + j_offset
                self.coordinates_per_node[node, 0] = grid_per_dim[0][i]
                self.coordinates_per_node[node, 1] = grid_per_dim[1][j]

        # -- set total number of cells and faces
        num_cells = nx * ny
        self.num_cells = num_cells
        self.num_faces = 4 * num_cells

        # -- constants for this mesh
        self.num_faces_per_cell = 4 * np.ones((num_cells), dtype=int)
        self.num_nodes_per_face = 2 * np.ones((4 * num_cells), dtype=int)

        # -- set nodes per face and faces per cell
        self.faces_per_cell = np.zeros((num_cells, 4), dtype=int)
        self.nodes_per_face = np.zeros((4 * num_cells, 2), dtype=int)
        for j in range(ny):
            for i in range(nx):
                # -- cell index
                cell = i + nx * j
                # -- faces per cell
                self.faces_per_cell[cell, 0] = 4 * cell
                self.faces_per_cell[cell, 1] = 4 * cell + 1
                self.faces_per_cell[cell, 2] = 4 * cell + 2
                self.faces_per_cell[cell, 3] = 4 * cell + 3
                # -- nodes per face (per cell)
                self.nodes_per_face[4 * cell, 0] = i + j * (nx + 1)
                self.nodes_per_face[4 * cell, 1] = i + 1 + j * (nx + 1)
                self.nodes_per_face[4 * cell + 1, 0] = i + 1 + j * (nx + 1)
                self.nodes_per_face[4 * cell + 1, 1] = i + 1 + (j + 1) * (nx + 1)
                self.nodes_per_face[4 * cell + 2, 0] = i + 1 + (j + 1) * (nx + 1)
                self.nodes_per_face[4 * cell + 2, 1] = i + (j + 1) * (nx + 1)
                self.nodes_per_face[4 * cell + 3, 0] = i + (j + 1) * (nx + 1)
                self.nodes_per_face[4 * cell + 3, 1] = i + j * (nx + 1)

        # -- enumerate boundary nodes per face ("side")
        # -- faces at x extrema
        nodes_per_side_xlow = np.zeros((ny, 2), dtype=int)
        nodes_per_side_xhig = np.zeros((ny, 2), dtype=int)
        for j in range(ny):
            nodes_per_side_xlow[ny - 1 - j, 0] = (j + 1) * (nx + 1)
            nodes_per_side_xlow[ny - 1 - j, 1] = j * (nx + 1)
            nodes_per_side_xhig[j, 0] = j * (nx + 1) + nx
            nodes_per_side_xhig[j, 1] = (j + 1) * (nx + 1) + nx
        # -- faces at y extrema
        nodes_per_side_ylow = np.zeros((nx, 2), dtype=int)
        nodes_per_side_yhig = np.zeros((nx, 2), dtype=int)
        for i in range(nx):
            nodes_per_side_ylow[i, 0] = i
            nodes_per_side_ylow[i, 1] = i + 1
            nodes_per_side_yhig[nx - 1 - i, 0] = i + 1 + ny * (nx + 1)
            nodes_per_side_yhig[nx - 1 - i, 1] = i + ny * (nx + 1)
        # -- compile into one side face array list (expeced as counter-clockwise in 2D)
        self.nodes_per_side = [nodes_per_side_ylow, nodes_per_side_xhig,
                               nodes_per_side_yhig, nodes_per_side_xlow]


# ------------------------------------------------------------------------------------------------ #
# orthogonal 3D mesh type
class orth_3d_mesh(base_mesh):
    '''
    Class for orthogonally structured 2D mesh data.
    This class generates an orthogonally structured mesh in an unstructured
    format suitable for creating unstructured-mesh input files.
    '''
    def __init__(self, bounds_per_dim, num_cells_per_dim):
        assert (False), 'orth_3d_mesh type not yet implemented...'


# ------------------------------------------------------------------------------------------------ #
# end of mesh_types.py
# ------------------------------------------------------------------------------------------------ #
