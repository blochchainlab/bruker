import numpy as np 

def main(inputfile, output):

	f = open(inputfile,'r')
	lines = f.readlines()
	f.close()

	vector_lines = [s for s in lines[21:-1]]
	vector_lines = [s[:-2] for s in vector_lines] # cut \n
	vector_str = [s.split(' ') for s in vector_lines]
	vector_list = [[float(i) for i in s] for s in vector_str]
	vector = np.array(vector_list)

	vector_norms = np.linalg.norm(vector, axis=1) 
	print('max = {}, >1 = {}'.format(vector_norms.max(), np.sum(vector_norms>1)))


	f2 = open(output,'w')

	for i in range(21):
		f2.writelines(lines[i])

	tot = len(vector_lines)
	for i in range(len(vector_lines)):
		v = vector[i]
		n = (np.round(v, 6)**2).sum()**0.5
		if n <= 1:
			f2.writelines('{:.6f} {:.6f} {:.6f}\n'.format(v[0], v[1], v[2]))
			tot -= 1
		else:
			v = v/n
			n = (np.round(v, 6)**2).sum()**0.5
			if n <= 1:
				f2.writelines('{:.6f} {:.6f} {:.6f}\n'.format(v[0], v[1], v[2]))
			else:
				v = v*(1-1e-6)
				n = (np.round(v, 6)**2).sum()**0.5
				if n <= 1:
					f2.writelines('{:.6f} {:.6f} {:.6f}\n'.format(v[0], v[1], v[2]))
				else:
					print('fucked')

	print('changed {} points'.format(tot))
	f2.writelines(lines[-1])
	f2.close()

if __name__ == "__main__":
	import sys
	main(*sys.argv[1:])
