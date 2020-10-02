#!/usr/bin/env python3
import sys
import numpy as np


def main(datalist, bvecfname, bvalfname):


	global_formater = '%.5f'
	print('\nUsing this output format: {}'.format(global_formater))


	## BVECS (PVM_DWGradVec + normalization)
	print('\nextracting bvecs from normalizing $PVM_DwGradVec')
	tag = '##$PVM_DwGradVec='
	bvec_data = grablines(datalist, tag)

	if len(bvec_data) >= 1:
		# read Number of dir
		Ndir = int(bvec_data[0].split('(')[1].split(',')[0].strip())

		# get dirs
		tmp = [float(i) for i in ' '.join(bvec_data[1:]).split(' ')]
		tmp = np.array(tmp).reshape((Ndir, 3))

		# attempt to detect b0
		b0_th = 1e-8
		b0_idx = np.where(np.all(np.abs(tmp) <= b0_th, axis=1))

		# normalize
		bvecs = tmp / np.linalg.norm(tmp, axis=1)[:,None]
		# fix div 0 error
		bvecs[np.isnan(bvecs)] = 0
		bvecs[np.isinf(bvecs)] = 0

		# save bvecs (in 3 columns, none of that 3 rows nonsense)
		print('Saving bvecs to {}'.format(bvecfname))
		np.savetxt(bvecfname, bvecs, fmt=global_formater)

	else:
		print('Couldn\'t create bvecs')


	## BVALS (PVM_DwEffBval)
	print('\nextracting bvals from $PVM_DwEffBval')
	tag = '##$PVM_DwEffBval='
	bval_data = grablines(datalist, tag)

	if len(bval_data) >= 1:
		# read Number of dir
		Ndir = int(bval_data[0].split('(')[1].split(')')[0].strip())

		# get bs
		tmp = [float(i) for i in ' '.join(bval_data[1:]).split(' ')]
		tmp = np.array(tmp)

		# warning for possibly weird b0s value
		if 'b0_idx' in locals():
			if len(b0_idx[0]) > 0:
				print('The detected b0s had the following bvalues:')
				print(tmp[b0_idx])

		# save bvals
		print('Saving bvals to {}'.format(bvalfname))
		np.savetxt(bvalfname, tmp, fmt=global_formater)

	else:
		print('Couldn\'t create bvals')



def findTagIndex(datalist, tag):
	# finds the first line of the desired ta in bruker method
	l = len(tag)
	idx = None
	for i in range(len(datalist)):
		if tag == datalist[i][:l]:
			idx = i
	return idx


def grablines(datalist, tag, newtag=['#','$']):
	# Grab all the lines
	# print from the tag to the the next newtag '#'
	idx = findTagIndex(datalist, tag)
	tmp = []
	if idx is None:
		print('Tag {} not found'.format(tag[3:-1]))
	else:
		for i in range(len(datalist[idx+1:])):
			flag = False
			for nt in newtag:
				if datalist[idx+1+i][0] == nt:
					flag = True
			if flag:
				idx2 = idx+1+i
				break
		for i in range(idx, idx2):
			tmp.append(datalist[i].rstrip())
	return tmp


if __name__ == "__main__":
	print('\nThis script create FSL style bvecs bvals from Bruker method')
	print('\nrequired parameters:')
	print('bruker_method_input bvec_output bval_output\n')

	# read method file
	f = open(sys.argv[1],'r')
	method_data = f.readlines()
	f.close()

	main(method_data, sys.argv[2], sys.argv[3])
