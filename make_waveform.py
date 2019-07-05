#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argparse
import numpy as np

from oct2py import Oct2Py

# make_waveform_param.m needs to be in the path
# make_waveform_param(gmax, tau, zeta, out_path, out_name, theta, phi, plotting, saving)
# defaults:
# gmax = 0.66 # mT/m, maximum gradient strenght [Magdeburg's gradient]
# tau = 10 # ms, gradient waveform duration, should be keep constant in an experiment even when changing the form, []
# zeta = np.math.acos(1/np.sqrt(3.)), half aperture of q cone, 0 ==> linear, 0 < zeta < acos(1/sqrt(3)) ==> cigar, acos(1/sqrt(3)) ==> sphere, acos(1/sqrt(3)) < zeta < pi/2 ==> pancake, pi/2 ==> planar
# out_name # output folder
# out_name # output name, Note that only ‘-‘ and ‘_’ special characters are allowed, extension in .gp
# theta, phi # Orientation of cone axis in lab frame
# plotting # boolean for displaying the waveform
# saving # boolean for saving the waveform file



DESCRIPTION = """
"""

EPILOG = """
"""

def buildArgsParser():
    p = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = DESCRIPTION,
    epilog = EPILOG)

    p._optionals.title = "Options and Parameters"

    p.add_argument(
        'gmax', action='store',
        metavar='gmax', type=float, default=0.66,
        help='maximum gradient strenght in T/m [%(default)s]')
    p.add_argument(
        'tau', action='store',
        metavar='tau', type=float, default=10e-3,
        help='waveform duration in s [%(default)s]')
    p.add_argument(
        'shape', action='store', metavar='shape', type=str,
        help='b-tensor shape, (lin, pro, sph, obl or sph)')
    p.add_argument(
        'outpath', action='store', metavar='outpath', type=str,
        help='waveform output folder path')
    p.add_argument(
        'outname', action='store', metavar='outname', type=str,
        help='waveform output name (ext .gp)')

    # this is ridiculous, it will be replaced by a dir file, just need to figure out which theta,phi convention this is
    p.add_argument(
        '--theta', dest='theta', action='store',
        metavar='theta', type=float, default=0.0,
        help='Orientation of cone axis in lab frame [%(default)s]')
    p.add_argument(
        '--phi', dest='phi', action='store',
        metavar='phi', type=float, default=0.0,
        help='Orientation of cone axis in lab frame [%(default)s]')
    p.add_argument(
        '-p', action='store_true', dest='plot', default=False,
        help='Plot waveform [%(default)s]')
    p.add_argument(
        '-s', action='store_false', dest='save', default=True,
        help='Don\'t save waveform [%(default)s]')

    return p




def main():

    parser = buildArgsParser()
    args = parser.parse_args()


    gmax = args.gmax
    if gmax >= 10:
        print('!!! maximum grad str (gmax) >= 10 T/m, this is suspicious, check your units')

    tau = args.tau
    if tau >= 100e-3:
        print('!!! waveform duration (tau) >= 100 ms, this is suspicious, check your units')
    elif tau <= 1e-3:
        print('!!! waveform duration (tau) <= 1 ms, this is suspicious, check your units')

    shape = args.shape
    # custom zeta eventually
    if shape == 'lin':
        zeta = 0
    elif shape == 'pro':
        zeta = np.acos(np.sqrt(2. / 3)) # Approximately 35 deg for prolate shape
    elif shape == 'sph':
        zeta = np.math.acos(np.sqrt(1 / 3.))
    elif shape == 'obl':
        zeta = 70. * np.pi/180. # Hardcoded to 70 deg, to achive oblate b-tensor
    elif shape == 'pla':
        zeta = np.pi/2.
    else:
        parser.error('btensor shape (shape) need to be one of lin, pro, sph, obl ,sph')
        return

    outpath = args.outpath
    outname = args.outname
    if outname[-3:] != '.gp':
        parser.error('output name (outname) need to end in .gp')
        return
    print('output = {}'.format(outpath+outname))

    theta = args.theta
    phi = args.phi
    print('orientation = ({},  {})'.format(theta, phi))
    
    plotting = args.plot
    if plotting:
        print('plot the waveform')
    saving = args.save
    if not saving:
        print('not saving the waveform')


    oc = Oct2Py()
    bb = oc.make_waveform(gmax, tau, zeta, outpath, outname, theta, phi, plotting, saving)

    # bvalue seems to scale with tau^3, which make sense because we KINd of have DELTA=delta
    # with normal trapeze pulse b ~ (G*delta)**2 * (DELTA-delta/3) ~ 2/3 G^2 delta^3
    print('b = {:.2e} s / m^2'.format(bb))
    print('b = {:.3f} ms / um^2'.format(bb*1e-9))

    # norms = np.linalg.norm(gxyz, axis=1)
    # print('maximum gradient norm = {}'.format(norms.max()))


if __name__ == "__main__":
    main()

