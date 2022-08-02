#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file   src/cdi_ipcress/python/plot_ipcress_opacity.py
# author Alex Long <along@lanl.gov>
# date   Monday, December 15, 2014, 5:44 pm
# brief  This script uses the functions in ipcress_reader.py to generate an interactive plot for
#        multigroup opacity data.
# note   Copyright (C) 2014-2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# import block
# ------------------------------------------------------------------------------------------------ #
import ipcress_reader as ip_reader
import matplotlib

from numpy import min, max
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from math import exp, log, floor, log10

import sys
if sys.version_info[0] < 3:
    from Tkinter import Tk, IntVar, DoubleVar, TOP, BOTH, Label, W, Entry, Radiobutton, Button, \
        BOTTOM, mainloop
else:
    from tkinter import Tk, IntVar, DoubleVar, TOP, BOTH, Label, W, Entry, Radiobutton, Button, \
        BOTTOM, mainloop
matplotlib.use('TkAgg')

# ------------------------------------------------------------------------------------------------ #
# main program
# ------------------------------------------------------------------------------------------------ #

# make sure an IPCRESS file is specified
if (len(sys.argv) != 2):
    print("Useage: {0} <path to ipcress file>".format(sys.argv[0]))
    sys.exit()

ipcress_file = sys.argv[1]

# set up TK plot
root = Tk()
root.wm_title("Data from IPCRESS File {0}".format(ipcress_file))

# font size for text in plots
# matplotlib.rcParams.update({'font.size': 18})

# get data dictionary from file
table_key_dict, materials = ip_reader.get_property_map_from_ipcress_file(ipcress_file)

# select the first material in the file as default and set as a TK variable
selected_ID = IntVar()  # make this a TK variable
selected_ID.set(materials[0])  # set TK variable

# add figuree
f = Figure(figsize=(5, 4), dpi=100)
a = f.add_subplot(111)

# set values for the first material (selected)
T_grid = table_key_dict["{0}_{1}".format("tgrid", selected_ID.get())]
rho_grid = table_key_dict["{0}_{1}".format("rgrid", selected_ID.get())]
hnu_grid = table_key_dict["{0}_{1}".format("hnugrid", selected_ID.get())]
mgr_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
mgp_grid = table_key_dict["{0}_{1}".format("pmg", selected_ID.get())]
mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]
mgrt_grid = table_key_dict["{0}_{1}".format("rtmg", selected_ID.get())]

# get log average rho and T for initial plot
avg_rho = exp(0.5 * (log(max(rho_grid)) + log(min(rho_grid))))
avg_rho = round(avg_rho, -int(floor(log10(avg_rho)) - 1))
avg_T = exp(0.5 * (log(max(T_grid)) + log(min(T_grid))))
avg_T = round(avg_T, -int(floor(log10(avg_T)) - 1))

# setup and initialize Tk values that can be tied to entry fields
# set temperature and density to log average value
target_rho = DoubleVar()
target_rho.set(avg_rho)
target_T = DoubleVar()
target_T.set(avg_T)
op_type = IntVar()
op_type.set(1)

# get interpolated opacity data for this temperature and density
# first check to see if the data is valid (non-zero)
mgr_valid = ip_reader.check_valid_data(mgr_grid)
mgp_valid = ip_reader.check_valid_data(mgp_grid)
mgs_valid = ip_reader.check_valid_data(mgs_grid)
mgrt_valid = ip_reader.check_valid_data(mgrt_grid)

name = selected_ID.get()
print("-------------------- BEGIN DATA PRINT FOR {0} ---------------------".format(name))
print("Group structure for {0} groups:".format(len(hnu_grid) - 1))
print(hnu_grid)
# if valid, interpolate data at target rho and target T
if (mgr_valid):
    mgr_interp = ip_reader.interpolate_mg_opacity_data(
        T_grid, rho_grid, hnu_grid, mgr_grid,
        target_rho.get(), target_T.get(), "multigroup absorption Rosseland")
if (mgp_valid):
    mgp_interp = ip_reader.interpolate_mg_opacity_data(
        T_grid, rho_grid, hnu_grid, mgp_grid,
        target_rho.get(), target_T.get(), "multigroup absorption Planckian")
if (mgs_valid):
    mgs_interp = ip_reader.interpolate_mg_opacity_data(
        T_grid, rho_grid, hnu_grid, mgs_grid,
        target_rho.get(), target_T.get(), "multigroup scattering")
if (mgrt_valid):
    # mgrt_interp =
    ip_reader.interpolate_mg_opacity_data(
        T_grid, rho_grid, hnu_grid, mgrt_grid,
        target_rho.get(), target_T.get(), "total rosseland")
print("-------------------- END DATA PRINT FOR {0} ---------------------".format(name))

# plotting data arrays
opr_data = []    # Rosseland absorption
opp_data = []    # Planck absorption
ops_data = []    # Rosseland scattering
oprt_data = []   # Rosseland total
hnu_data = []    # hnu plot data

# add points to plotting data twice to get the bar structure and find max and min opacity
min_opacity = 1.0e10
max_opacity = 0.0
for hnu_i, hnu in enumerate(hnu_grid[:-1]):
    hnu_data.append(hnu)
    hnu_data.append(hnu_grid[hnu_i + 1])
    if (mgr_valid):
        opr_data.append(mgr_interp[hnu_i])
        opr_data.append(mgr_interp[hnu_i])
        min_opacity = min([min(opr_data), min_opacity])
        max_opacity = max([max(opr_data), max_opacity])

    if (mgp_valid):
        opp_data.append(mgp_interp[hnu_i])
        opp_data.append(mgp_interp[hnu_i])
        min_opacity = min([min(opp_data), min_opacity])
        max_opacity = max([max(opp_data), max_opacity])

    if (mgs_valid):
        ops_data.append(mgs_interp[hnu_i])
        ops_data.append(mgs_interp[hnu_i])
        min_opacity = min([min(ops_data), min_opacity])
        max_opacity = max([max(ops_data), max_opacity])

# Plot data with the initial settings
a.set_xscale('log')
a.set_yscale('log')
a.set_xlabel('hnu (keV)')
a.set_ylabel('opacity (cm sq./g)')

if (mgr_valid):
    a.plot(hnu_data, opr_data, 'b-', label="{0} Rosseland Absorption".format(name))
else:
    print("ERROR: Invalid multigroup Rosseland absorption data (all zeros)")

if (mgs_valid):
    a.plot(hnu_data, ops_data, 'r-', label="{0} Rosseland Scattering".format(name))
else:
    print("ERROR: Invalid multigroup scattering data (all zeros)")

a.set_xlim([0.9 * min(hnu_grid), 1.1 * max(hnu_grid)])
a.set_ylim([0.7 * min_opacity, 1.3 * max_opacity])

a.legend(loc='best')

canvas = FigureCanvasTkAgg(f, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


# ------------------------------------------------------------------------------------------------ #
# Plot loop executed when the "Plot" button is pressed
# ------------------------------------------------------------------------------------------------ #
def plot_op():
    f.clf()

    a = f.add_subplot(111)

    # make sure selected material is available, if not found set to first material
    mat_found = False
    for mat in materials:
        if selected_ID.get() == mat:
            mat_found = True
    if mat_found is False:
        selected_ID.set(materials[0])

    # get values associated with this ID from dictionary
    T_grid = table_key_dict["{0}_{1}".format("tgrid", selected_ID.get())]
    rho_grid = table_key_dict["{0}_{1}".format("rgrid", selected_ID.get())]
    hnu_grid = table_key_dict["{0}_{1}".format("hnugrid", selected_ID.get())]
    mgr_grid = table_key_dict["{0}_{1}".format("ramg", selected_ID.get())]
    mgp_grid = table_key_dict["{0}_{1}".format("pmg", selected_ID.get())]
    # scattering is only availabe with Rosseland weighting
    mgs_grid = table_key_dict["{0}_{1}".format("rsmg", selected_ID.get())]
    mgrt_grid = table_key_dict["{0}_{1}".format("rtmg", selected_ID.get())]

    # get interpolated opacity data for this temperature and density
    # first check to see if the data is valid (non-zero)
    mgr_valid = ip_reader.check_valid_data(mgr_grid)
    mgp_valid = ip_reader.check_valid_data(mgp_grid)
    mgs_valid = ip_reader.check_valid_data(mgs_grid)
    mgrt_valid = ip_reader.check_valid_data(mgrt_grid)

    name = selected_ID.get()
    print("-------------------- BEGIN DATA PRINT FOR {0} --------------------".format(name))
    # if valid, interpolate data at target rho and target T
    if (mgr_valid):
        mgr_interp = ip_reader.interpolate_mg_opacity_data(
            T_grid, rho_grid, hnu_grid, mgr_grid, target_rho.get(), target_T.get(),
            "multigroup absorption Rosseland")
    if (mgp_valid):
        mgp_interp = ip_reader.interpolate_mg_opacity_data(
            T_grid, rho_grid, hnu_grid, mgp_grid, target_rho.get(), target_T.get(),
            "multigroup absorption Planckian")
    if (mgs_valid):
        mgs_interp = ip_reader.interpolate_mg_opacity_data(
            T_grid, rho_grid, hnu_grid, mgs_grid, target_rho.get(), target_T.get(),
            "multigroup scattering")
    if (mgrt_valid):
        # mgrt_interp
        ip_reader.interpolate_mg_opacity_data(
            T_grid, rho_grid, hnu_grid, mgrt_grid, target_rho.get(), target_T.get(),
            "total rosseland")
    print("-------------------- END DATA PRINT FOR {0} --------------------".format(name))

    # plotting data arrays
    opr_data = []    # Rosseland absorption
    opp_data = []    # Planck absorption
    ops_data = []    # Rosseland scattering
    # oprt_data = []   # Rosseland total
    hnu_data = []    # hnu plot data

    # add points to plotting data twice to get the bar structure and find max and min opacity
    min_opacity = 1.0e10
    max_opacity = 0.0
    for hnu_i, hnu in enumerate(hnu_grid[:-1]):
        hnu_data.append(hnu)
        hnu_data.append(hnu_grid[hnu_i + 1])

        if (mgr_valid):
            opr_data.append(mgr_interp[hnu_i])
            opr_data.append(mgr_interp[hnu_i])
            min_opacity = min([min(opr_data), min_opacity])
            max_opacity = max([max(opr_data), max_opacity])

        if (mgp_valid):
            opp_data.append(mgp_interp[hnu_i])
            opp_data.append(mgp_interp[hnu_i])
            min_opacity = min([min(opp_data), min_opacity])
            max_opacity = max([max(opp_data), max_opacity])

        if (mgs_valid):
            ops_data.append(mgs_interp[hnu_i])
            ops_data.append(mgs_interp[hnu_i])
            min_opacity = min([min(ops_data), min_opacity])
            max_opacity = max([max(ops_data), max_opacity])

    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlabel('hnu (keV)')
    a.set_ylabel('opacity (cm sq./g)')

    # use label for Planck or Rosseland
    if (op_type.get() == 1 and mgr_valid):
        a.plot(hnu_data, opr_data, 'b-', label="{0} Rosseland Absorption".format(name))
    elif (op_type.get() == 2 and mgp_valid):
        a.plot(hnu_data, opp_data, 'g-', label="{0} Planckian Absorption".format(name))
    elif (op_type.get() == 3 and mgr_valid and mgp_valid):
        a.plot(hnu_data, opr_data, 'b-', label="{0} Rosseland Absorption".format(name))
        a.plot(hnu_data, opp_data, 'g-', label="{0} Planckian Absorption".format(name))
    else:
        if (not mgr_valid):
            print("ERROR: Invalid multigroup Rosseland absorption data (all zeros)")
        if (not mgp_valid):
            print("ERROR: Invalid multigroup Planckian absorption data (all zeros)")
    if (mgs_valid):
        a.plot(hnu_data, ops_data, 'r-', label="{0} Rosseland Scattering".format(name))
    else:
        print("ERROR: Invalid multigroup scattering data (all zeros)")

    a.legend(loc='best')
    a.set_xlim([0.9 * min(hnu_grid), 1.1 * max(hnu_grid)])
    a.set_ylim([0.7 * min_opacity, 1.3 * max_opacity])
    canvas.draw()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)


# ------------------------------------------------------------------------------------------------ #
# set up entry fields and their associated labels
# ------------------------------------------------------------------------------------------------ #
L_rho = Label(root, text="Density ({0} g/cc - {1} g/cc)".format(min(rho_grid), max(rho_grid)))
L_rho.pack(anchor=W)
target_rho_field = Entry(root, textvariable=target_rho)
target_rho_field.pack(anchor=W)

L_T = Label(root, text="Temperature ({0} keV - {1} keV)".format(min(T_grid), max(T_grid)))
L_T.pack(anchor=W)
target_T_field = Entry(root, textvariable=target_T)
target_T_field.pack(anchor=W)

material_str = ""
for mat in materials:
    material_str = "{0} {1}".format(material_str, mat)
L_material = Label(root, text="Materials in file: {0} ".format(material_str))
L_material.pack(anchor=W)
material_field = Entry(root, textvariable=selected_ID)
material_field.pack(anchor=W)


# ------------------------------------------------------------------------------------------------ #
def _quit():
    root.quit()     # stops mainloop
    # Avoid python error on Win32: 'Fatal Python Error: PyEval_RestoreThread: NULL tstate'
    root.destroy()


# ------------------------------------------------------------------------------------------------ #
L_opacity = Label(root, text="Opacity Weighting Method")
L_opacity.pack(anchor=W)
Radiobutton(root, text="Rosseland", padx=20, variable=op_type, value=1).pack(anchor=W)
Radiobutton(root, text="Planckian",    padx=20, variable=op_type, value=2).pack(anchor=W)
Radiobutton(root, text="Show Both",    padx=20, variable=op_type, value=3).pack(anchor=W)

button = Button(master=root, text='Quit', command=_quit)
button.pack(side=BOTTOM)

button = Button(root, text='Plot', command=plot_op)
button.pack(side=BOTTOM)

mainloop()

# ------------------------------------------------------------------------------------------------ #
# End plot_ipcress_opacity.py
# ------------------------------------------------------------------------------------------------ #
