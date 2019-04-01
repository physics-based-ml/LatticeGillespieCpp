import os

# Directory idxs
i_start = 1 # inclusive
i_end = 300 # inclusive

for i in range(i_start,i_end+1):
	dir_name = "lattice_v%03d" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%03d/counts/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%03d/lattice/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%03d/nns/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)