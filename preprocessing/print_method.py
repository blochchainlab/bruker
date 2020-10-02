#!/usr/bin/env python3
import sys
import numpy as np


def findTagIndex(datalist, tag):
	l = len(tag)
	idx = None
	for i in range(len(datalist)):
		if tag == datalist[i][:l]:
			idx = i
	return idx

def print1(datalist, tag):
	# print the tag's line
	idx = findTagIndex(datalist, tag)
	if idx is None:
		print('No {} found\n'.format(tag[3:-1]))
	else:
		print(datalist[idx])

def print2(datalist, tag):
	# print the tag's line and the next
	idx = findTagIndex(datalist, tag)
	if idx is None:
		print('No {} found\n'.format(tag[3:-1]))
	else:
		print(datalist[idx].rstrip())
		print(datalist[idx+1])

def printmany(datalist, tag, newtag=['#','$']):
	# print from the tag to the the next newtag '#'
	idx = findTagIndex(datalist, tag)
	if idx is None:
		print('No {} found\n'.format(tag[3:-1]))
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
			print(datalist[i].rstrip())
		print('\n')


def printmethod(datalist):
	# Print an handpicked list of Bruker tag and values from the method file

	tag1 = ['##$Method=',
			'##$PVM_RareFactor=',
			'##$IROffset=',
			'##$IRSpacing=',
			'##$NIRPoints=',
			'##$PVM_EchoTime=',
			'##$PVM_RepetitionTime=',
			'##$PVM_NAverages=',
			'##$PVM_ScanTime=',
			'##$PVM_SliceThick=',
			'##$PVM_NRepetitions=',
			'##$PVM_DwBvalEach=',
			'##$PVM_DwAoImages=',
			'##$PVM_DwEffBval=',
			'##$PVM_DwNDiffDir=',
			'##$PVM_DwNDiffExpEach=',
			'##$PVM_DiffPrepMode=']


	tag2 = ['##$PVM_ScanTimeStr',
			'##$EffectiveTE=',
			'##$PVM_FovCm=',
			'##$PVM_SpatResol=',
			'##$PVM_SPackArrNSlices=',
			'##$PVM_Matrix=',
			'##$PVM_DwGradDur=',
			'##$PVM_DwGradSep']

	tagmany = ['##$PVM_SliceGeo=',
			   '##$PVM_DwDir=',
			   '##$PVM_DwGradVec=',
			   '##$PVM_DwEffBval=']


	for t in tag1:
		print1(datalist, t)
	for t in tag2:
		print2(datalist, t)
	for t in tagmany:
		printmany(datalist, t)


	print('\n\nHeuristic Voxel check\n')
	fovcm = np.array([float(i) for i in datalist[findTagIndex(datalist, '##$PVM_FovCm=')+1].rstrip().split(' ')])
	resmm = np.array([float(i) for i in datalist[findTagIndex(datalist, '##$PVM_SpatResol=')+1].rstrip().split(' ')])
	matrix = np.array([float(i) for i in datalist[findTagIndex(datalist, '##$PVM_Matrix=')+1].rstrip().split(' ')])

	print('resmm = 10*fovcm / matrix = {}\n'.format(10*fovcm / matrix))
	print('resmm={} fovcm={} matrix={}'.format(resmm,fovcm,matrix))




if __name__ == "__main__":
	
	# read method file
	f = open(sys.argv[1],'r')
	method_data = f.readlines()
	f.close()

	printmethod(method_data)

