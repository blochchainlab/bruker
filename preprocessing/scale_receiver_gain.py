#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import argparse

import nibabel as nib
import numpy as np


desciption = """
Loads a Nifti file for receiver gain intensity correction
The RG can either be fetched from the acqp file produced by the Bruker scan or directly supplied
The outputed data will have intensity original/RG and will be saved in float32
"""

# TODO:
# check that the pattern finding didnt fail


def _build_args_parser():
    p = argparse.ArgumentParser(description=desciption, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('data', metavar='data', help='Path of the nifti file.')
    # need either a pah to acqp or a number for receiver gain
    groupRG = p.add_mutually_exclusive_group(required=True)
    groupRG.add_argument('--acqp')
    groupRG.add_argument('--val')
    p.add_argument('output', metavar='output', help='Path of the output nifti.')
    return p

def main():
    parser = _build_args_parser()
    args = parser.parse_args()

    img_data = nib.load(args.data)
    affine = img_data.affine


    # acqp has priority
    if args.acqp is None:
        RG = float(args.val)
    else:
        # read file
        f= open(args.acqp, 'r')
        l = f.readlines()
        f.close()
        # find the line
        pattern = "##$RG="
        idx = np.where([ll[:len(pattern)] == pattern for ll in l])[0][0]
        # read value
        RG = float(l[idx][len(pattern):].strip())


    print('Receiver gain is {}'.format(RG))

    # prepare scaled data
    scaled_data = img_data.get_fdata() / RG

    # save as float32 at output
    nib.nifti1.Nifti1Image(scaled_data.astype(np.float32), affine).to_filename(args.output)

    print('data divided by RG save at {}'.format(args.output))


if __name__ == "__main__":
    main()

