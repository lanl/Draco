#!/usr/bin/env python3
# -------------------------------------------*-python-*------------------------------------------- #
# file   src/predict/python/convert_pt_to_nn.py
# author Mathew Cleveland
# date   February 22nd 2022
# brief  This script converts py-torch *.pt files to a simple draco Neural Network format
#
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

import torch
from torchinfo import summary
import struct
import argparse


def simple_write(model, filename, binary=False):
    info = summary(model, verbose=2, col_names=["kernel_size"])
    n_layers = 0
    for layer in info.summary_list[1:]:
        if ('weight' in layer.inner_layers):
            n_layers += 1

    if binary:
        Activations = {'None': 0, 'ReLU': 1}
        Layer_Type = {'Linear': 0}
        file = open(filename, 'wb')
        file.write(bytearray(struct.pack("i", 1)))
        file.write(bytearray(struct.pack("i", n_layers)))
        last_activation = "None"
        for layer in info.summary_list[1:]:
            if ('weight' not in layer.inner_layers):
                last_activation = layer.class_name
            else:
                file.write(bytearray(struct.pack("i", Activations[last_activation])))
                file.write(bytearray(struct.pack("i", Layer_Type[layer.class_name])))
                shape = layer.inner_layers['weight']['kernel_size'].strip('[]').split(",")
                file.write(bytearray(struct.pack("i", int(shape[0]))))
                file.write(bytearray(struct.pack("i", int(shape[1]))))
                shape = layer.inner_layers['bias']['kernel_size'].strip('[]')
                file.write(bytearray(struct.pack("i", int(shape))))
                last_activation = "None"
        # write out the parameters
        for p in model.parameters():
            if p.requires_grad:
                for vals in p.data:
                    try:
                        for val in vals:
                            file.write(bytearray(struct.pack("f", float(val))))
                    except Exception:
                        file.write(bytearray(struct.pack("f", float(vals))))
    else:
        file = open(filename, 'w')
        # Write out the layout
        file.write("Layers: " + str(n_layers) + "\n")
        last_activation = "None"
        for layer in info.summary_list[1:]:
            if ('weight' not in layer.inner_layers):
                last_activation = layer.class_name
            else:
                file.write(last_activation + " " + layer.class_name + " "
                           + layer.inner_layers['weight']['kernel_size'] + " "
                           + layer.inner_layers['bias']['kernel_size'] + "\n")
                last_activation = "None"
        # write out the parameters
        for p in model.parameters():
            if p.requires_grad:
                file.write('\n')
                for vals in p.data:
                    file.write(str(vals.tolist()) + '\n')


def parse_arguments():
    """Read arguments from a command line."""
    parser = argparse.ArgumentParser(
        description='Convert py - torch NN * .pt file to a simple draco NN file * .nnb'
        ' ( or *.nn for non-binary text format)')
    parser.add_argument('-pt', dest="pt_files", required=True, nargs="+",
                        help='Torch files to convert to binary format')
    args = parser.parse_args()
    return args


def main(args):
    args = parse_arguments()
    for pt_file in args.pt_files:
        net = torch.jit.load(pt_file)
        name = pt_file.replace('.pt', '.nnb')
        print("Writing binary file -> " + name)
        simple_write(net, name, True)
        name = pt_file.replace('.pt', '.nn')
        print("Writing ASCII file -> " + name)
        simple_write(net, name, False)
    pass


if __name__ == '__main__':
    args = parse_arguments()
    main(args)
