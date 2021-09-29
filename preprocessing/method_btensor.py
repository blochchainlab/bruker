#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import argparse
# import pylab as pl # TODO
from gnlc_waveform.btensor import compute_q_from_G, compute_B_from_q, get_btensor_eigen, get_btensor_shape_topgaard
# from gnlc_waveform.viz import peraxisplot, plot # TODO





DESCRIPTION = """
Parse method files to create btensor list and associated files
Intended for the Budde Bruker sequence
"""

EPILOG = """
Michael Paquette
29.09.2021
"""

def buildArgsParser():
    p = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = DESCRIPTION,
    epilog = EPILOG)

    # p._optionals.title = "Options and Parameters"

    p.add_argument(
        'fname', action='store', metavar='method_fname', type=str,
        help='Path to the method file')

    p.add_argument(
        'out', action='store', metavar='output_basename', type=str,
        help='Basename for output')

    # TODO
    # p.add_argument(
    #     '-p', action='store_true', dest='plot', default=False,
    #     help='Plot waveform [%(default)s]')

    return p





# read all lines and parse into list
def rawtext(fname):
    with open(fname, 'r') as f:
        return f.readlines()


# find a tag and grap everything until the next '$$' or '##'
def grab_tag(rawmethod, tag):
    # remove tag marker if there
    if tag[:3] == '##$':
        tag = tag[3:]
    # loop list until found
    tagindex = None
    for i,l in enumerate(rawmethod):
        if l[:len(tag)+3] == '##$' + tag:
            tagindex = i
            break
    if i is None:
        print('tag not found')
        return None
    # find the end of tag
    end_tag = None
    j = 1
    while end_tag is None:
        tmp = rawmethod[i+j][:2]
        if (tmp == '$$') or (tmp == '##'):
            end_tag = i+j
        j += 1
    
    return rawmethod[i:end_tag]


# parse single line float tag
def parse_tag_single_line_float(rawtag):
    return np.array([float(rawtag[0].strip().split('=')[-1])])


# parse multiline float tag
def parse_tag_multi_line_float(rawtag):
    arrayshape = tuple([int(n) for n in rawtag[0].strip().split('=')[-1][1:-1].split(',')])
    array1D = np.concatenate([np.array([float(n) for n in l.strip().split(' ')]) for l in rawtag[1:]])
    return array1D.reshape(arrayshape)


# parse multiline string tag
def parse_tag_multi_line_string(rawtag):
    return rawtag[1].strip()[1:-1]


# parse any tag using heuristic to decide between:
#     parse_tag_single_line_float
#     parse_tag_multi_line_float
#     parse_tag_multi_line_string
def parse_tag(rawtag):
    if len(rawtag) == 1:
        return parse_tag_single_line_float(rawtag)
    elif rawtag[1][0] == '<':
        return parse_tag_multi_line_string(rawtag)
    else:
        return parse_tag_multi_line_float(rawtag)
    

# TODO handle waveform where A and B have different dt
# waveforms to btensor
def W2B(waveA, waveB, tA, tpause, G):
    # waveform A and waveform B should have norm <= 1
    # waveform B is inverted, as if 180degree pulse
    # tA is waveform A lenght in s
    # waveform B assumes to be same lenght and discretization
    # tpause is he waveform separation in s
    # G is the maximum gradient strenght in T/m
    dt = tA / (waveA.shape[0]-1) # s
    npause = int(np.floor(tpause / dt))
    gt = np.concatenate((waveA, np.zeros((npause, 3)), -waveB), axis=0)
    gradient = gt*G
    qt = compute_q_from_G(gradient, dt)
    btensor = compute_B_from_q(qt, dt)
    eigval, eigvec = get_btensor_eigen(btensor)
    bval, bs, bp, bl = get_btensor_shape_topgaard(eigval)

    #        s/m^2,   s/m^2,     N/A, s/m^2, N/A, N/A, N/A 
    return btensor,  eigval,  eigvec,  bval,  bs,  bp,  bl







def main():

    parser = buildArgsParser()
    args = parser.parse_args()

    method_fname = args.fname
    baseout = args.out

    # load a method file of interest
    rawmethod = rawtext(method_fname)


    # tags of interest
    tags = ['EchoTime',              # Sequence param (ms)
            'PVM_RepetitionTime',    # Sequence param (ms)
            'DwGradDur1',            # Waveform A duration (ms)
            'DwGradDur2',            # Waveform B duration (ms)
            'DwGradTsep',            # Between A and B duration (ms)
            'DwDir',                 # bvecs
            'DwDynGradShapeEnum1',   # Name waveform A
            'DwDynGradShapeEnum2',   # Name waveform B
            'DwGradAmpG',            # Gradient amplitude of each bval (mT/m)
            'DwGradAmpScale',        # Gradient amplitude of each bval (percent of max)      
            'DwGradShapeArray1',     # Waveform A
            'DwGradShapeArray2',     # Waveform B      
            'DwGradAmpScalarArray',  # Gradient amplitude of each volume (percent of max)
            'DwR00',                 # Element (0,0) of rotation matrix
            'DwR11',                 # Element (1,1) of rotation matrix        
            'DwR22',                 # Element (2,2) of rotation matrix
            'DwR01',                 # Element (0,1) of rotation matrix
            'DwR02',                 # Element (0,2) of rotation matrix
            'DwR10',                 # Element (1,0) of rotation matrix
            'DwR12',                 # Element (1,2) of rotation matrix
            'DwR20',                 # Element (2,0) of rotation matrix
            'DwR21',                 # Element (2,1) of rotation matrix
            'DwDirsFilename',        # Name of bvec file
            ]


    # grab all tags of interest
    rawtags = {}
    for tag in tags:
        rawtags[tag] =  grab_tag(rawmethod, tag)


    # parse all tags of interest
    parsedtag = {}
    for k in rawtags.keys():
        parsedtag[k] = parse_tag(rawtags[k])


    # Print sequence parameters
    print('\nSequence Summary')
    print('TR = {:.2f} ms'.format(parsedtag['PVM_RepetitionTime'][0]))
    print('TE = {:.2f} ms'.format(parsedtag['EchoTime'][0]))
    
    print('\nWaveA filename = {:}'.format(parsedtag['DwDynGradShapeEnum1']))
    print('WaveA = {:.2f} ms'.format(parsedtag['DwGradDur1'][0]))
    print('WaveB filename = {:}'.format(parsedtag['DwDynGradShapeEnum2']))
    print('WaveB = {:.2f} ms'.format(parsedtag['DwGradDur2'][0]))
    print('Wave Separation = {:.2f} ms'.format(parsedtag['DwGradTsep'][0]))

    print('\nbvec filename = {:}'.format(parsedtag['DwDirsFilename']))
    print('Nb bvec = {:.0f}'.format(parsedtag['DwDir'].shape[0]))

    print('\nNb bval = {:.0f}'.format(parsedtag['DwGradAmpG'].shape[0]))
    print('G = {:} mT/m'.format(['{:.2f}'.format(f) for f in parsedtag['DwGradAmpG']]))
    print('G = {:} %'.format(['{:.1f}'.format(f) for f in parsedtag['DwGradAmpScale']]))

    print('\nNb volume = {:.0f}'.format(parsedtag['DwGradAmpScalarArray'].shape[0]))


    # build rotation matrix
    # these matrix bring [0,0,1] to each bvec
    rot_mats = np.zeros((parsedtag['DwR00'].shape[0], 3, 3))
    rot_mats[:, 0, 0] = parsedtag['DwR00']
    rot_mats[:, 0, 1] = parsedtag['DwR01']
    rot_mats[:, 0, 2] = parsedtag['DwR02']
    rot_mats[:, 1, 0] = parsedtag['DwR10']
    rot_mats[:, 1, 1] = parsedtag['DwR11']
    rot_mats[:, 1, 2] = parsedtag['DwR12']
    rot_mats[:, 2, 0] = parsedtag['DwR20']
    rot_mats[:, 2, 1] = parsedtag['DwR21']
    rot_mats[:, 2, 2] = parsedtag['DwR22']


    Gmax = parsedtag['DwGradAmpG'][0] / (0.01*parsedtag['DwGradAmpScale'][0])
    print('\nScanner Gmax = {:.2f} mT/m'.format(Gmax))

    # unrotated waveform with Gmax
    btensor_init, eigval_init, eigvec_init, bval_init, bs_init, bp_init, bl_init = W2B(parsedtag['DwGradShapeArray1'][0], 
                                                                                       parsedtag['DwGradShapeArray2'][0], 
                                                                                       parsedtag['DwGradDur1']*1e-3, 
                                                                                       parsedtag['DwGradTsep']*1e-3, 
                                                                                       Gmax*1e-3)


    print('initial waveform btensor with Gmax (s/mm^2)')
    print(btensor_init*1e-6)
    print(' ')

    print('eigen val and vec')
    print('{:.1f} s/mm^2'.format(eigval_init[0]*1e-6))
    print('[{:.3f} {:.3f} {:.3f}]'.format(*eigvec_init[:,0]))
    print(' ')
    print('{:.1f} s/mm^2'.format(eigval_init[1]*1e-6))
    print('[{:.3f} {:.3f} {:.3f}]'.format(*eigvec_init[:,1]))
    print(' ')
    print('{:.1f} s/mm^2'.format(eigval_init[2]*1e-6))
    print('[{:.3f} {:.3f} {:.3f}]'.format(*eigvec_init[:,2]))
    print(' ')

    print('Normalized shape (Spherical, Planar, Linear)')
    print('[{:.2f} {:.2f} {:.2f}]'.format(bs_init, bp_init, bl_init))


    btensors = np.zeros((rot_mats.shape[0], 3, 3))
    eigenvec_max = np.zeros((rot_mats.shape[0], 3))
    eigenvec_min = np.zeros((rot_mats.shape[0], 3))
    eigenvals = np.zeros((rot_mats.shape[0], 3))
    bvals = np.zeros(rot_mats.shape[0])
    normalized_shape = np.zeros((rot_mats.shape[0], 3))

    for i in range(rot_mats.shape[0]):
        rotmat = rot_mats[i]
        G = 0.01*parsedtag['DwGradAmpScalarArray'][i] * Gmax

        # rotate waveform
        btensor_rot, eigval_rot, eigvec_rot, bval_rot, bs_rot, bp_rot, bl_rot = W2B(parsedtag['DwGradShapeArray1'][0].dot(rotmat.T), 
                                                                                    parsedtag['DwGradShapeArray2'][0].dot(rotmat.T), 
                                                                                    parsedtag['DwGradDur1']*1e-3, 
                                                                                    parsedtag['DwGradTsep']*1e-3, 
                                                                                    G*1e-3)

        # clean
        eigval_rot = np.clip(eigval_rot, 0, np.inf)
        bs_rot = np.clip(bs_rot, 0, 1)
        bp_rot = np.clip(bp_rot, 0, 1)
        bl_rot = np.clip(bl_rot, 0, 1)


        # print('[{:.3f} {:.3f} {:.3f}]'.format(*eigval_rot*1e-6)) # eigenval
        # print('[{:.3f} {:.3f} {:.3f}]'.format(*eigvec_rot[:,0])) # eigenvec of least eigenval

        # print(btensor_rot*1e-6)
        btensors[i] = btensor_rot*1e-6 # s/mm^2

        # print(eigvec_rot[:, np.argmax(eigval_rot)])
        # print(np.linalg.norm(eigvec_rot[:, np.argmax(eigval_rot)]))
        eigenvec_max[i] = eigvec_rot[:, np.argmax(eigval_rot)]

        # print(eigvec_rot[:, np.argmin(eigval_rot)])
        # print(np.linalg.norm(eigvec_rot[:, np.argmin(eigval_rot)]))
        eigenvec_min[i] = eigvec_rot[:, np.argmin(eigval_rot)]

        # print(bval_rot*1e-6)
        bvals[i] = bval_rot*1e-6 # s/mm^2

        # print(eigval_rot*1e-6)
        eigenvals[i] = eigval_rot*1e-6 # s/mm^2

        # print(bs_rot, bp_rot, bl_rot)
        normalized_shape[i] = bs_rot, bp_rot, bl_rot


    print('\nSaving flattened btensors (in s/mm^2) as...')
    print(baseout + '_btensor.txt')
    np.savetxt(baseout + '_btensor.txt', btensors.reshape((-1, 9)))

    print('\nSaving biggest eigenvec as...')
    print(baseout + '_e1.txt')
    np.savetxt(baseout + '_e1.txt', eigenvec_max)

    print('\nSaving smallest eigenvec as...')
    print(baseout + '_e3.txt')
    np.savetxt(baseout + '_e3.txt', eigenvec_min)

    print('\nSaving eigenvals (in s/mm^2) as...')
    print(baseout + '_eigval.txt')
    np.savetxt(baseout + '_eigval.txt', eigenvals)

    print('\nSaving bvals (in s/mm^2) as...')
    print(baseout + '_bvals.txt')
    np.savetxt(baseout + '_bvals.txt', bvals)

    print('\nSaving normalized Westin shape (Sph, Pla, Lin) as...')
    print(baseout + '_wshape.txt')
    np.savetxt(baseout + '_wshape.txt', normalized_shape)



if __name__ == "__main__":
    main()


