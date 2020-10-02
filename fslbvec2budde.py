import numpy as np 


def main(inputfile, outfolder, outputfile):
	bvec = np.genfromtxt(inputfile)
	if bvec.shape[1] != 3:
		bvec = bvec.T
	# bvec = bvec.T

	f = open(outfolder+outputfile, 'w')

	f.writelines('# {}\n'.format(outputfile))

	f.writelines('#This file contains a table of diffusion gradient vectors\n')
	f.writelines('# to be loaded for custom methods.\n')
	f.writelines('# Each set should have a single number on a line\n')
	f.writelines('# followed by the subesequent number of directions\n')
	f.writelines('# Any number of directions is allowed.\n')
	f.writelines('# But, no vectors with a norm > |1| are allowed.\n')
	f.writelines('#\n')
	f.writelines('#\n')
	f.writelines('#\n')

	f.writelines('{}\n'.format(bvec.shape[0]))

	for i in range(bvec.shape[0]):
		x,y,z = bvec[i]
		f.writelines('{:.6f} {:.6f} {:.6f}\n'.format(x, y, z))

	f.close()


if __name__ == "__main__":
	print('usage: fsl.bvec outfolder output.dat')
	print('assumes single shell')
	import sys
	main(*sys.argv[1:])

