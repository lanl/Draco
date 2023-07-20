#!/usr/bin/env python
# -------------------------------------------*-python-*------------------------------------------- #
# file   src/cdi_ipcress/python/ipcress_reader.py
# author Alex Long <along@lanl.gov>
# date   Monday, December 15, 2014, 5:44 pm
# brief This script has functions that parse an IPCRESS file and returns a dictionary that contains
#        data for each property and each material present in the file. This script also contains
#        interpolation functions for opacity data.
# note   Copyright (C) 2016-2023 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# import block
# ------------------------------------------------------------------------------------------------ #
import numpy as np
from struct import pack, unpack
from math import log, exp

# ------------------------------------------------------------------------------------------------ #
# These are the functions that are used to read data from the binary IPCRESS file. It also contains
# a function for interpolating in density and temperature. The data locations are specified in
# cdi_ipcress/doc/IPCRESS_File_Format.pdf
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def get_data_for_id(filename, data_start_index, num_entries):
    temp_grid = []
    # "rb" is read binary mode
    with open(filename, "rb") as f:
        f.seek(data_start_index * 8)
        for i in range(num_entries):
            word = f.read(8)
            temp_grid.append(unpack('>d', word)[0])
    return temp_grid


# ------------------------------------------------------------------------------------------------ #
def write_data_for_id(filename, data_start_index, num_entries, new_values):
    # "wb" is write binary mode
    with open(filename, "r+b") as f:
        f.seek(data_start_index * 8)
        for i in range(num_entries):
            s = pack('>d', new_values[i])
            f.write(s)


# ------------------------------------------------------------------------------------------------ #
def interpolate_opacity_date_checks(target_rho, target_T, rho_grid, T_grid, print_str=""):
    # don't allow extrapolation
    if (target_rho < np.min(rho_grid)):
        target_rho = np.min(rho_grid)
    if (target_rho > np.max(rho_grid)):
        target_rho = np.max(rho_grid)
    if (target_T < np.min(T_grid)):
        target_T = np.min(T_grid)
    if (target_T > np.max(T_grid)):
        target_T = np.max(T_grid)
    if (print_str is not None):
        print("Interpolating {0}--Target rho: {1} , target T: {2}".format(
            print_str, target_rho, target_T))


# ------------------------------------------------------------------------------------------------ #
def interpolate_mg_opacity_data(T_grid, rho_grid, hnu_grid, op_data, target_rho, target_T,
                                print_str=""):
    n_rho = len(rho_grid)
    n_hnu = len(hnu_grid)

    # don't allow extrapolation
    interpolate_opacity_date_checks(target_rho, target_T, rho_grid, T_grid, print_str)

    # get correct index of adjacent density points
    rho_L = 1000
    rho_G = 0
    for rho_i, rho in enumerate(rho_grid[:-1]):
        if (target_rho >= rho and target_rho <= rho_grid[rho_i + 1]):
            rho_L = rho_i
            rho_G = rho_i + 1
            break

    # get correct index of adjacent temperature points
    T_L = 1000
    T_G = 0
    for T_i, T in enumerate(T_grid[:-1]):
        if (target_T >= T and target_T <= T_grid[T_i + 1]):
            T_L = T_i
            T_G = T_i + 1
            break

    # get the adjacent rows of the opacity index
    # get the points of the opacity index
    rho_L_T_L = op_data[n_rho * T_L * (n_hnu - 1) + rho_L * (n_hnu - 1):n_rho * T_L * (n_hnu - 1) +
                        rho_L * (n_hnu - 1) + (n_hnu - 1)]
    rho_L_T_G = op_data[n_rho * T_G * (n_hnu - 1) + rho_L * (n_hnu - 1):n_rho * T_G * (n_hnu - 1) +
                        rho_L * (n_hnu - 1) + (n_hnu - 1)]
    rho_G_T_L = op_data[n_rho * T_L * (n_hnu - 1) + rho_G * (n_hnu - 1):n_rho * T_L * (n_hnu - 1) +
                        rho_G * (n_hnu - 1) + (n_hnu - 1)]
    rho_G_T_G = op_data[n_rho * T_G * (n_hnu - 1) + rho_G * (n_hnu - 1):n_rho * T_G * (n_hnu - 1) +
                        rho_G * (n_hnu - 1) + (n_hnu - 1)]

    interp_op = []
    # interpolate for each frequency point
    for i in range(n_hnu - 1):
        # if (rho_L_T_L[i] < 1.0e-10) or (rho_L_T_G[i] < 1.0e-10) or (rho_G_T_L[i] < 1.0e-10) or
        # (rho_G_T_G[i] < 1.0e10):
        #    interp_op.append(1.0e-10)
        # print("{0} {1} {2} {3}" .format(rho_L_T_L[i], rho_L_T_G[i], rho_G_T_L[i], rho_G_T_G[i]))
        log_op_T_L = log(rho_L_T_L[i]) + log(target_rho / rho_grid[rho_L]) / \
            log(rho_grid[rho_G] / rho_grid[rho_L]) * log(rho_G_T_L[i] / rho_L_T_L[i])
        log_op_T_G = log(rho_L_T_G[i]) + log(target_rho / rho_grid[rho_L]) / \
            log(rho_grid[rho_G] / rho_grid[rho_L]) * log(rho_G_T_G[i] / rho_L_T_G[i])
        log_op = log_op_T_L + log(target_T / T_grid[T_L]) / log(T_grid[T_G] / T_grid[T_L]) * \
            (log_op_T_G - log_op_T_L)
        interp_op.append(exp(log_op))

    if (print_str is not None):
        print("hnu(keV)            opacity(sq_cm/g)         opacity(1/cm)")
        for i, hnu in enumerate(hnu_grid[:-1]):
            print("{0}     {1}     {2}".format(0.5 * (hnu + hnu_grid[i + 1]), interp_op[i],
                                               interp_op[i] * target_rho))
    return interp_op


# ------------------------------------------------------------------------------------------------ #
def interpolate_gray_opacity_data(T_grid, rho_grid, op_data, target_rho, target_T, print_str=""):
    n_rho = len(rho_grid)

    # don't allow extrapolation
    interpolate_opacity_date_checks(target_rho, target_T, rho_grid, T_grid, print_str)

    rho_L = 1000
    rho_G = 0
    for rho_i, rho in enumerate(rho_grid[:-1]):
        if (target_rho >= rho and target_rho <= rho_grid[rho_i + 1]):
            rho_L = rho_i
            rho_G = rho_i + 1
            break

    for T_i, T in enumerate(T_grid[:-1]):
        if (target_T >= T and target_T <= T_grid[T_i + 1]):
            T_L = T_i
            T_G = T_i + 1
            break

    # get the adjacent rows of the opacity index
    rho_L_T_L = op_data[n_rho * T_L + rho_L]
    rho_L_T_G = op_data[n_rho * T_G + rho_L]
    rho_G_T_L = op_data[n_rho * T_L + rho_G]
    rho_G_T_G = op_data[n_rho * T_G + rho_G]

    # interpolate in log space
    log_op_T_L = log(rho_L_T_L) + log(target_rho / rho_grid[rho_L]) / \
        log(rho_grid[rho_G] / rho_grid[rho_L]) * log(rho_G_T_L / rho_L_T_L)
    log_op_T_G = log(rho_L_T_G) + log(target_rho / rho_grid[rho_L]) / \
        log(rho_grid[rho_G] / rho_grid[rho_L]) * log(rho_G_T_G / rho_L_T_G)
    log_op = log_op_T_L + log(target_T / T_grid[T_L]) / \
        log(T_grid[T_G] / T_grid[T_L]) * (log_op_T_G - log_op_T_L)
    interp_op = exp(log_op)

    return interp_op


# ------------------------------------------------------------------------------------------------ #
def read_information_from_file(ipcress_file):

    word_array = []
    with open(ipcress_file, "rb") as f:
        for i in range(26):
            word = f.read(8)
            if not word:
                break
            word_array.append(word)
            # print(int(unpack('>d', word)[0]))

    # title = word_array[0]
    toc_int = []
    offset = 2
    for i in range(offset, offset + 24):
        toc_int.append(int(unpack('>d', word_array[i])[0]))

    n_data_records = toc_int[14]
    # mxrec = toc_int[1] - toc_int[0]
    mxkey = toc_int[16]
    # print("Number of data records: {0}".format(n_data_records))
    # print("Beginnging of data: {0}".format(toc_int[0]))
    # print("Max records: {0} , max search keys: {1}".format(mxrec, mxkey))

    mat_property = []

    ds = []
    dfo = []
    # tdf = []
    num_mats = 0
    mat_ids = []
    with open(ipcress_file, "rb") as f:
        # Read in array that lists the data sizes in this file
        f.seek(toc_int[0] * 8)
        # print("Table of data sizes")
        for i in range(n_data_records):
            word = f.read(8)
            ds.append(int(unpack('>d', word)[0]))

        # Read in array gives the offsets between data
        f.seek(toc_int[1] * 8)
        # print("Table of data file offesets")
        for i in range(n_data_records):
            word = f.read(8)
            dfo.append(int(unpack('>d', word)[0]))

        # Read in material IDs present in this file
        f.seek(dfo[0] * 8)
        # print("Table of    material identifiers")
        word = f.read(8)
        num_mats = int(unpack('>d', word)[0])
        for i in range(num_mats):
            word = f.read(8)
            mat_ids.append(int(unpack('>d', word)[0]))

        # Read in list of properties in this file available for each material
        # entries in this table are 24 bytes each
        f.seek(toc_int[10] * 8)
        # print("Table of data fields for each    material")
        word = f.read(72)  # ignore the first 72 bytes, they don't contain useful information
        for i in range(1, toc_int[14]):
            # oredering is "matID" "data type" "fill"
            temp_property = []
            for j in range(mxkey):
                three_string = []
                three_string.append(f.read(8).decode("utf-8"))
                three_string.append(f.read(8).decode("utf-8"))
                three_string.append(f.read(8).decode("utf-8"))
                if (j == 0):
                    temp_property.append(three_string[2].strip())
                elif (j == 1):
                    temp_property.append(three_string[0].strip())
                else:
                    temp_property.append(i)  # index of data table containing values
            try:
                temp_property = [temp_property[0].decode('ascii'),
                                 temp_property[1].decode('ascii'), temp_property[2]]
                mat_property.append(temp_property)
            except:  # noqa: E722
                mat_property.append(temp_property)

    materials = []
    for m in range(num_mats):
        materials.append([m, mat_ids[m]])

    # return the list of available properties, data file offsets and data sizes
    return materials, mat_property, dfo, ds


# ------------------------------------------------------------------------------------------------ #
def write_information_to_file(ipcress_file, material_ID, mat_property, new_values):
    materials, property_list, dfo, ds = read_information_from_file(ipcress_file)
    # check to make sure material is in file
    material_IDs = []
    for imat in materials:
        material_IDs.append(str(imat[1]))
    if (not (material_ID in material_IDs)):
        print("ERROR: Material ID not found in file, not changing anything!")
        return

    # try to find property in file
    property_found = False
    # propery_index = 0
    for prop_i, prop in enumerate(property_list):
        if (material_ID == prop[0] and mat_property == prop[1]):
            property_found = True
            property_index = prop_i
            break

    # make sure sizes match of property you're about to write
    if (property_found and ds[property_index + 1] != len(new_values)):
        print("ERROR: Number of new values does not match size of old values, not changing "
              "anything!")
        return

    # if the combination of property and material was found, write the new data to the ipcress file
    if property_found:
        write_data_for_id(ipcress_file, dfo[property_index + 1], ds[property_index + 1], new_values)
    else:
        print("ERROR: Combination of material ID and property not found, not changing anything!")
    return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Checks to see if there are any zeros in the opcaity data--zero data is difficult to handle and for
# now we are going to ignore data sets that contain zeros and print an error message
def check_valid_data(opacity_grid):
    for item in opacity_grid:
        if (item != 0.0):
            return True
    return False
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# return a dictionary where the keys are "<material ID>_<property_ID>" and the values are the data
def get_property_map_from_ipcress_file(ipcress_file):
    # load data from IPCRESS file
    # dfo is the array of data file offsets, ds is the array of data sizes
    materials, property_list, dfo, ds = read_information_from_file(ipcress_file)

    # build dictionary of data, keys are "property_matID"
    table_key_dict = {}
    for prop_i, prop in enumerate(property_list):
        table_key_dict["{0}_{1}".format(prop[1], prop[0])] = get_data_for_id(
            ipcress_file, dfo[prop_i + 1], ds[prop_i + 1])

    material_list = []
    for material in materials:
        material_list.append(material[1])

    return table_key_dict, material_list

# ------------------------------------------------------------------------------------------------ #
# End ipcress_reader.py
# ------------------------------------------------------------------------------------------------ #
