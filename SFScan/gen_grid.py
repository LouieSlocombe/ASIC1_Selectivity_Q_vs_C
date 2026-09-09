#!/usr/bin/env python3

# Define arguments

import sys, getopt

directory = ''
system = ''
cutoff=1.5
radius=5.5
figure=''
output = ''
try:
   opts, args = getopt.getopt(sys.argv[1:],"hd:s:c:r:f:o:",["dir=","sys=","cut=","rad=","--fig=","out="])
except getopt.GetoptError as err:
   print(err)
   print('ERROR: Syntaxe should be "gen_grid.py -d <directory> -s <coord_file> -c <cutoff> -r <radius> -f <figure> -o <output>"')
   sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
      print('gen_grid.py -d <directory> -s <coord_file> PDB -c <cutoff> -r <radius> -f <figure> PNG -o <output>')
      sys.exit()
   elif opt in ("-d", "--dir"):
      directory = arg
   elif opt in ("-s", "--sys"):
      system = arg
   elif opt in ("-c", "--cut"):
      cutoff = arg
   elif opt in ("-r", "--rad"):
      radius = arg
   elif opt in ("-f", "--fig"):
      figure = arg
   elif opt in ("-o", "--out"):
      output = arg

##########################

import numpy as np
import plotly.graph_objects as go
import molparse as mp
import time
import os
import mout

sf_sys = mp.parse(f'{directory}/{system}')
sf_sys.summary()

HSD = sf_sys['rHSD']
GLY = sf_sys['rGLY']
ALA = sf_sys['rALA']
SER = sf_sys['rSER']

whole_sys_com = np.array(sf_sys.CoM()[:2])
# print(whole_sys_com)

coms = []

for res in HSD:

	# res.summary()
	print(res.index,res.number,res.CoM(verbosity=0))
	coms.append(res.CoM(verbosity=0))

for res in GLY:

	# res.summary()
	print(res.index,res.number,res.CoM(verbosity=0))
	coms.append(res.CoM(verbosity=0))

for res in ALA:

	# res.summary()
	print(res.index,res.number,res.CoM(verbosity=0))
	coms.append(res.CoM(verbosity=0))

for res in SER:

	# res.summary()
	print(res.index,res.number,res.CoM(verbosity=0))
	coms.append(res.CoM(verbosity=0))

max_x = max([v[0] for v in coms])
min_x = min([v[0] for v in coms])
max_y = max([v[1] for v in coms])
min_y = min([v[1] for v in coms])

print(f'{max_x=}')
print(f'{min_x=}')
print(f'{max_y=}')
print(f'{min_y=}')

grid_points = []

ION = sf_sys.atoms[-1]
# print(ION.position)
z_coordinate = ION.np_pos[-1]
# print(z_coordinate)

cutoff=float(cutoff)
radius=float(radius)
spacing = 0.5

x_grid = np.arange(min_x,max_x,spacing)
# print(x_grid)
y_grid = np.arange(min_y,max_y,spacing)
# print(y_grid)

for x in x_grid:
	for y in y_grid:

		if np.linalg.norm(whole_sys_com - np.array([x,y])) >= radius:
			continue

		nearby = False
		for atom in sf_sys.atoms[:-2]:
			
			if np.linalg.norm(atom.np_pos - np.array([x,y,z_coordinate])) <= cutoff:
				nearby = True
				break

		if nearby:
			continue

		grid_points.append([x,y,z_coordinate])

# print(grid_points[0:3])
n_points = len(grid_points)
print(f'{n_points=}')

fig = sf_sys.plot3d(flat=False,show=False)

trace = go.Scatter3d(x=[v[0] for v in grid_points],y=[v[1] for v in grid_points],z=[v[2] for v in grid_points],mode='markers')

fig.add_trace(trace)

# fig.update_yaxes(
#     scaleanchor = "x",
#     scaleratio = 1
# )

mp.write(f'{directory}/{figure}',fig)

# Generate PDB strucutre for each points in the grid
for i, coordinate in enumerate(grid_points):

	start = time.perf_counter()

	output_name = f'{output}_{i:03}.pdb'

	ION.position = grid_points[i]

	mp.write(f'{directory}/{output_name}',sf_sys,verbosity=0)

	end = time.perf_counter()

	mout.progress(i,len(grid_points),prepend=f'{i=:03}',append=f' {end-start:.1f} sec')

mout.progress(i,i,prepend=f'{i=:03}',append='         DONE')