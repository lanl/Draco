# Mesh File Generation Python Subdirectory (src/mesh/python)

This subdirectory, `src/mesh/python/`, is intended to serve as a location for
different methods of generating X3D or RTT mesh files.
The X3D and RTT mesh file formats are supported by C++ parsers in
`src/RTT_Format_Reader/` and `src/mesh/`
There is a python module called `mesh_types` with mesh classes, where each
class has its own `__init__` method for calculating data needed to create X3D
files (some of which should be usable for RTT file creation).
For instance, one could add a mesh class that stochastically samples points and
triangulates them to form a mesh of triangles.
Currently, there is a script called `x3d_generator.py`, which takes command-line
input to instantiate a mesh object from one of the mesh_type classes and outputs
a set of X3D files (main mesh file and boundary node files) of the mesh object.

## Example Usage

To create X3D files for a 4x4x4 3D orthogonal structured mesh on domain (x,y,z)
in [0,1]x[0,2]x[0,4]:

```bash
./x3d_generator.py --mesh_type orth_3d_mesh --num_per_dim 4 4 4 --bnd_per_dim 0 1 0 2 0 4
```
