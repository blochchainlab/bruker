#!/usr/bin/env python3
import numpy as np
import nibabel as nib
import pylab as pl
import argparse



desciption = """
Correct the mean instensity of the input volume to be equal inside a brain mask.
Uses the minimum as the normalization target

"""




def _build_args_parser():
    p = argparse.ArgumentParser(description=desciption, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('data', metavar='data', help='Path of the nifti file.')
    p.add_argument('output', metavar='output', help='Path of the output nifti.')
    p.add_argument('scale', metavar='scale', help='Path of the output scaling.')
    p.add_argument('--mask', metavar='mask', help='Path of the brain mask for normalization.')
    return p


def main():
    parser = _build_args_parser()
    args = parser.parse_args()

    # data
    img_data = nib.load(args.data)
    data = img_data.get_fdata()
    affine = img_data.affine

    # mask
    mask_path = args.mask
    if mask_path is None:
        print('No mask given, using the full volumes')

        # mean intensity inside each volume
        volume_mean_intensity = data.mean(axis=(0,1,2))
    else:
        img_mask = nib.load(args.mask)
        mask = img_mask.get_fdata().astype(np.bool)

        # mean intensity inside mask for each volume
        volume_mean_intensity = data[mask].mean(axis=0)

    # desired mean inensity of each volume
    # we use the min as a proxy for steady-state
    target_intensity = volume_mean_intensity.min()


    pl.figure()
    pl.plot(volume_mean_intensity, color='blue', label='mean intensity in mask')
    pl.axhline(target_intensity, color='red', alpha=0.25, label='target for renomalization')
    pl.legend()
    pl.show()


    # since we use min() as target, scalings are <= 1
    scalings = target_intensity/volume_mean_intensity

    # scaled data
    scaled_data = data.copy() * scalings[None,None,None,:]

    # save as float32 at output
    nib.nifti1.Nifti1Image(scaled_data.astype(np.float32), affine).to_filename(args.output)
    print('Scaled data saved at {}'.format(args.output))

    # save scalings
    np.savetxt(args.scale, scalings)
    print('Scalings saved at {}'.format(args.scale))



if __name__ == "__main__":
    main()






