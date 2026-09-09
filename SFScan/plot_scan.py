#!/usr/bin/env python3

'''
	
	To-Do's
	=======

	- Determine size of binding pocket by considering the area of the heatmap below a certain threshold
	- RMSD between two transformed heatmaps

	Analysis
	========

	- Take the 6 sodium replicas where you substitute the ion

		> report the change in binding area
		> report the symmetry/asymmetry
		> look at the location of the minima
		> some description of the depth of the binding pocket
		> look at the closest interacting atoms from the minimum (compare to coordination number)

'''

scan_start = 0.9
scan_end = 5.0
scan_step = 0.1
interpolation_steps = 200

recalculate = False

import sys
import os
import glob
import plotly.graph_objects as go
import numpy as np
import mout
from scipy.interpolate import griddata, CloughTocher2DInterpolator, LinearNDInterpolator, interp1d
from molparse import hijack
import molparse as mp
from ase.optimize import BFGS
from ase.optimize import MDMin
from ase.neb import NEB, interpolate
import pickle

import plotly.io as pio   
pio.kaleido.scope.mathjax = None

mout.hideDebug()

class Scan(object):

	def __init__(self,key,rst_key=None,wall_height=10,locut=0.9,hicut=2.0,path='.',cap=None,mid_walls=False):
		self.key = key
		self.rst_key = rst_key
		self.wall_height = wall_height
		self.hicut = hicut
		self.locut = locut
		self.path = path
		self.cap = cap

		self.mid_walls = mid_walls

		self.reference_position = None

		self.missing = []
		self.notconverged = []
		self.samples = []
		self.pes_slices = []

		self._grid = None
		self._transformed_points = None
		self.neb_images = None

		self.transformation_matrix = None

		self.decay_path_rc1 = None
		self.decay_path_rc2 = None
		self.decay_path_rc3 = None
		self.decay_path_rc4 = None

		self.example_sys = None

		self.minimum = None
		self.min_e = None

		self.read_energy_and_pdbs()

		if rst_key is not None:
			self.read_separation_json()

	def get_reference_position(self):

		rep = self.key.split('/')
		file = f'{self.key}/sf{rep[-1][3]}.pdb'

		try:
			sys = mp.parse(file)
		except:
			mout.warningOut(f'No {file} for scan {self.key}')
			return

		self.reference_position = sys.atoms[-1].position[:2]

	def read_separation_json(self):
		import json
		self.sep_dict = json.load(open(f"{self.rst_key}/sep.json"))
		print(f"{self.rst_key}/sep.json")
		self.bb_sep = self.sep_dict['DGN9-DCN1']
		print(self.sep_dict['DGN9-DCN1'])

	def read_energy_and_pdbs(self):

		pdb_files = sorted(glob.glob(f'{self.path}/{self.key}/scan_???.pdb'))
		mout.varOut("Number of pdb files",len(pdb_files))

		energies = []
		positions = []
		missing = []

		for i,pdb in enumerate(pdb_files):
			mout.progress(i,len(pdb_files))

			energy_files = glob.glob(pdb.replace('.pdb','.energy'))

			# if not empty list
			if energy_files:

				# get the QM region energy
				with open(energy_files[0]) as f:
					for line in f:
						energy = float(line.split()[-1]) # in Hartrees
						energies.append(energy*630)
						break

				# get the coordinate of the ion
				sys = mp.parse(pdb,verbosity=0)
				ion = sys.atoms[-1]
				positions.append(ion.position[:2])

				self.samples.append(ion.position[:2])

				if self.example_sys is None:
					self.example_sys = sys
					self.get_histidine_coms()
					self.get_glycine_coms()
					# print('\r',end='')
					self.example_sys.remove_chain('I',verbosity=0)

			else:
				# print('\r',end='')
				mout.warningOut(f"{pdb} has no associated energy file!")
				# get the coordinate of the ion
				sys = mp.parse(pdb,verbosity=0)
				ion = sys.atoms[-1]
				missing.append(ion.position[:2])

				self.notconverged.append(ion.position[:2])

		mout.progress(len(pdb_files),len(pdb_files))

		self.locut = min([min(p) for p in positions])
		self.hicut = max([max(p) for p in positions])

		print(f'{self.locut=}')
		print(f'{self.hicut=}')
		print(f'Missing {len(missing)} energy files.')

		# find the energy minimum
		self.min_e = min(energies)

		# shift the energies so they are relative to the minimum
		energies = [e - self.min_e for e in energies]
		
		self.energy_cutoff = max(energies)

		# reshape the data
		data = [[x,y,e] for (x,y),e in zip(positions,energies)]

		self.point_list = data

	def plot_heatmap(self,show=False,show_samples=False,contours=True,transformed=False,show_missing=True):

		# we have a list of [x,y,e] points
		if transformed:
			interpolator = self.get_interpolator(True)

			nx, ny = (100, 100)
			xv = np.linspace(self.transformed_locut, self.transformed_hicut, nx)
			yv = np.linspace(self.transformed_locut, self.transformed_hicut, ny)

		else:

			interpolator = self.get_interpolator(False)

			nx, ny = (100, 100)
			xv = np.linspace(self.locut, self.hicut, nx)
			yv = np.linspace(self.locut, self.hicut, ny)

		mesh_data = []
		for y in yv:

			y_slice = []
			for x in xv:
				energy = interpolator(x,y)

				if energy > self.energy_cutoff:
					y_slice.append(None)
				else:
					y_slice.append(energy)

			mesh_data.append(y_slice)

		mesh_data = np.array(mesh_data)

		fig = go.Figure()

		if show_samples:
			if transformed:
				trace = go.Scatter(name='samples',x=[s[0] for s in self.transform(self.samples)],y=[s[1] for s in self.transform(self.samples)],mode='markers')
				fig.add_trace(trace)
			else:
				trace = go.Scatter(name='samples',x=[s[0] for s in self.samples],y=[s[1] for s in self.samples],mode='markers')
				fig.add_trace(trace)

		if show_missing:

			p = self.notconverged
			if transformed:
				p = self.transform(p)
			x = [m[0] for m in p]
			y = [m[1] for m in p]

			trace = go.Scatter(name='not converged',x=x,y=y,mode='markers',marker_symbol='x',marker_color='red')
			fig.add_trace(trace)

		if self.pes_slices:
			for i,(p,z) in enumerate(self.pes_slices):
				if transformed:
					p = self.transform(p)
				x = [v[0] for v in p]
				y = [v[1] for v in p]
				trace = go.Scatter(name=f'slice {i}',x=x,y=y,mode='lines')
				fig.add_trace(trace)

		if self.reference_position is not None:
			x,y = self.reference_position
			if transformed:
				x,y = self.transform([x,y])			
			trace = go.Scatter(name='ion reference',x=[x],y=[y],mode='markers')
			fig.add_trace(trace)

		if contours:
			contours = dict(start=0.0,end=self.energy_cutoff,size=self.energy_cutoff/20)
			trace = go.Contour(name=self.key,x=xv,y=yv,z=mesh_data,contours=contours,showscale=False)
		else:
			trace = go.Heatmap(name=self.key,x=xv,y=yv,z=mesh_data)
		
		fig.add_trace(trace)

		if not transformed:
			fig = self.example_sys.plot3d(fig=fig,flat=True,show=False)
		else:
			fig = self.example_sys.plot3d(fig=fig,flat=True,show=False,transform=self.transform)

		if self.minimum:
			x=self.minimum.rcs[0]
			y=self.minimum.rcs[1]

			if transformed:
				x,y = self.transform([x,y])

			trace = go.Scatter(name='minimum',x=[x],y=[y],mode='markers')
			fig.add_trace(trace)

		fig.update_layout(title=self.key)

		if show:
			fig.show()

		return fig

	def extend_walls(self):
		mout.out(f'placing walls {self.wall_height}, {self.locut - 0.2}:{self.hicut + 0.2}')
		self.point_list.append([self.locut - 0.2,self.locut - 0.2,self.wall_height])
		self.point_list.append([self.locut - 0.2,self.hicut + 0.2,self.wall_height])
		self.point_list.append([self.hicut + 0.2,self.hicut + 0.2,self.wall_height])
		self.point_list.append([self.hicut + 0.2,self.locut - 0.2,self.wall_height])
		
		if self.mid_walls:
			self.point_list.append([self.locut - 0.2,(self.locut + self.hicut)/2,self.wall_height])
			self.point_list.append([self.hicut + 0.2,(self.locut + self.hicut)/2,self.wall_height])
			self.point_list.append([(self.locut + self.hicut)/2,self.locut - 0.2,self.wall_height])
			self.point_list.append([(self.locut + self.hicut)/2,self.hicut + 0.2,self.wall_height])

	def transform(self,point):
		if isinstance(point[0], list):
			result = []
			for p in point:
				result.append(self.transform(p))
			return result
		assert len(point) == 2
		if point[0] is None and point[1] is None:
			return [None,None]
		else:
			u,v,_ = np.matmul(self.transformation_matrix,np.array([*point,1]).T)
		return [u,v]

	@property
	def transformed_points(self):
		if self._transformed_points is None:

			# closest histidine
			h_index = np.argmin([np.linalg.norm(np.array(h) - np.array(self.minimum.rcs)) for h in self.histidine_coms])
			other_indices = [0,1,2]
			other_indices.pop(h_index)

			# transformation matrix (https://stackoverflow.com/questions/55546892/transform-points-from-one-triangle-to-another-triangle)
			A = np.array([
							[0,np.sqrt(3)/2,np.sqrt(3)/2],
							[0,0.5,-0.5],
							[1,1,1],
					     ])

			B = np.array([
							[self.histidine_coms[h_index][0],self.histidine_coms[other_indices[0]][0],self.histidine_coms[other_indices[1]][0]],
							[self.histidine_coms[h_index][1],self.histidine_coms[other_indices[0]][1],self.histidine_coms[other_indices[1]][1]],
							[1,1,1],
					     ])

			self.transformation_matrix = np.matmul(A,np.linalg.inv(B))

			self._transformed_points = []
			for x,y,e in self.point_list:
				u,v = self.transform([x,y])
				self._transformed_points.append([u,v,e])

		self.transformed_locut = min([min(p[:2]) for p in self._transformed_points])
		self.transformed_hicut = max([max(p[:2]) for p in self._transformed_points])

		return self._transformed_points
	
	def get_interpolator(self,transformed=False,absolute=False):
		shift = 0.0
		if absolute:
			shift = self.min_e
		if transformed:
			self.interpolator = CloughTocher2DInterpolator([[p[0], p[1]] for p in self.transformed_points],[p[2]+shift for p in self.point_list])
		else:
			self.interpolator = CloughTocher2DInterpolator([[p[0], p[1]] for p in self.point_list],[p[2]+shift for p in self.point_list])
		return self.interpolator

	def get_calculator(self):
		return hijack.FakeSurfaceCalculator(pmf=self.get_interpolator())

	def optimize(self,atoms):

		atoms.set_calculator(self.get_calculator())

		dyn = BFGS(atoms=atoms)#,logfile=f'{self.key}_{atoms.name}_opt.log')

		mout.headerOut(f"Running optimiser ({atoms.name})...")
		dyn.run()

	def neb(self,n_images,k,reactant,product,fmax=0.04):

		self.neb_images = [reactant.copy() for i in range(n_images//2)] + [product.copy() for i in range(n_images//2)]

		interpolate(self.neb_images)

		neb = NEB(self.neb_images,k=k)

		for image in self.neb_images[1:-1]:
			# image.set_constraints(locut=0.91,hicut=2.28)
			image.set_calculator(self.get_calculator())

		optimizer = BFGS(neb, trajectory=f'{self.key}_neb.traj', logfile=f'{self.key}_neb.log')

		mout.headerOut(f"Running NEB... logfile={self.key}_neb.log")

		try:
			optimizer.run(fmax=fmax,steps=2000)
		except hijack.NaNEncounteredInPmfError:
			mout.errorOut("NEB did not finish correctly!")

		import scipy.signal as sps
		self.neb_minima_indices = [0] + list(sps.argrelextrema(np.array(self.get_neb_energies()), np.less)[0]) + [-1]
		self.neb_maxima_indices = list(sps.argrelextrema(np.array(self.get_neb_energies()), np.greater)[0])

		self.neb_stationary_x = [self.get_neb_rcdist()[i] for i in self.neb_minima_indices + self.neb_maxima_indices]
		self.neb_stationary_y = [self.get_neb_energies()[i] for i in self.neb_minima_indices + self.neb_maxima_indices]

		self.neb_stationary_points = [[x,y] for x,y in zip(self.get_neb_rcdist(),self.get_neb_energies())]

	def get_neb_energies(self):
		interpolator = self.get_interpolator()
		return [interpolator(x,y) for x,y in zip(*self.get_neb_rcpath())]

	def get_neb_rcdist(self,normalise=None):

		# renormalise to distance between canonical and tautomer

		path = [[x,y] for x,y in zip(*self.get_neb_rcpath())]

		distances = [0.0]
		for i,p in enumerate(path[1:]):
			d = np.array(p) - np.array(path[i])
			distances.append(distances[-1]+np.linalg.norm(d))

		if normalise is not None:
			distances = list(normalise*np.array(distances)/distances[-1])

		return distances

	def get_neb_rcpath(self,walls=None):

		if walls is not None:

			images = self.neb_images

			# canonical

			can_wall_images = []

			vec = np.array(images[0].rcs) - np.array(images[1].rcs)

			for i in range(walls[0]):
				rcs = list(np.array(images[0].rcs) + (i+1)*vec)
				can_wall_images.append(rcs)

			# tautomer

			tau_wall_images = []

			vec = np.array(images[-1].rcs) - np.array(images[-2].rcs)

			for i in range(walls[1]):
				rcs = list(np.array(images[-1].rcs) + (i+1)*vec)
				tau_wall_images.append(rcs)

			# print(reversed(can_wall_images),[i.rcs for i in self.neb_images],tau_wall_images)

			combined = list(reversed(can_wall_images)) + [i.rcs for i in self.neb_images] + tau_wall_images

			return [p[0] for p in combined], [p[1] for p in combined]

		else:
			return [i.rcs[0] for i in self.neb_images], [i.rcs[1] for i in self.neb_images]

	def get_gmx_like_neb_path(self,walls=None):

		gmx_rc = []

		x,y = self.get_neb_rcpath(walls=walls)

		for rc1,rc2 in zip(x,y):
			gmx_rc.append((rc1+rc2)/2)

		return gmx_rc

	def plot_energy_slices(self,show=False):
		
		fig = go.Figure()

		self.pes_slices = []

		# find energies along linear paths between histidine CoM's

		path = [self.histidine_coms[0],self.histidine_coms[1]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='HSD0 --> HSD1',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.H0_H1_slice_min = min(z)

		path = [self.histidine_coms[0],self.histidine_coms[2]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='HSD0 --> HSD2',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.H0_H2_slice_min = min(z)

		path = [self.histidine_coms[1],self.histidine_coms[2]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='HSD1 --> HSD2',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.H1_H2_slice_min = min(z)

		path = [self.glycine_coms[0],self.glycine_coms[1]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='GLY0 --> GLY1',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.G0_G1_slice_min = min(z)

		path = [self.glycine_coms[0],self.glycine_coms[2]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='GLY0 --> GLY2',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.G0_G2_slice_min = min(z)

		path = [self.glycine_coms[1],self.glycine_coms[2]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='GLY1 --> GLY2',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		self.G1_G2_slice_min = min(z)

		path = [self.minimum.rcs,self.histidine_coms[0]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='min --> HSD0',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		path = [self.minimum.rcs,self.histidine_coms[1]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='min --> HSD1',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		path = [self.minimum.rcs,self.histidine_coms[2]]
		p, z = self.get_pes_slice(path,steps=50)
		fig.add_trace(go.Scatter(name='min --> HSD2',x=[i-len(z)/2 for i,_ in enumerate(z)],y=z))
		
		fig.add_hline(self.interpolator(*self.minimum.rcs))
		# fig.add_trace(go.Scatter(name='minimum',x=[25],y=[self.interpolator(*self.minimum.rcs)]))
		
		if show: 
			fig.show()

		mp.write(f'{self.key}/slices.pdf',fig)

	def get_histidine_coms(self):

		self.histidine_coms = []

		for residue in self.example_sys.residues['HSD']:

			nitrogen = residue.atoms['N'][2]

			self.histidine_coms.append(nitrogen.position[:2])

		return self.histidine_coms

	def get_glycine_coms(self):

		self.glycine_coms = []

		for residue in self.example_sys.residues['GLY']:

			# residue.summary()
			if residue.atoms['O']:
				oxygen = residue.atoms['O'][0]

				self.glycine_coms.append(oxygen.position[:2])

		return self.glycine_coms

	def get_pes_slice(self,path,steps=None):

		xdata = []
		ydata = []
		zdata = []

		if steps:

			new_path = []

			for i,(x,y) in enumerate(path):

				if i == 0:
					continue

				new_path.append(path[i-1])

				for j in range(steps):

					this_x = (path[i][0] - path[i-1][0]) * (j+1)/steps + path[i-1][0]
					this_y = (path[i][1] - path[i-1][1]) * (j+1)/steps + path[i-1][1]

					new_path.append([this_x,this_y])

			path = new_path

		for x,y in path:

			xdata.append(x)
			ydata.append(y)

			e = float(self.interpolator(x,y))
			
			if np.isnan(e):
				continue

			zdata.append(e)

		self.pes_slices.append([path,zdata])

		return path, zdata

def main():

	if recalculate:

		keys = [
			# 'Na/REP0', # ok
			# 'Na/REP1', # ok
			# 'Na/REP2', # ok
			# 'Na/REP3', # ok
			# 'Na/REP4', # ok
			# 'Na/REP5', # ok

			# 'Na_HW/REP0', # ok
			# 'Na_HW/REP1', # ok
			# 'Na_HW/REP2', # ok
			# 'Na_HW/REP3', # ok
			# 'Na_HW/REP4', # ok
			# 'Na_HW/REP5', # ok

			# 'K/REP0', # ok
			# 'K/REP1', # ok
			# 'K/REP2', # ok
			# 'K/REP3', # ok
			# 'K/REP4', # ok
			# 'K/REP5', # ok

			# 'K_HW/REP0', # ok
			# 'K_HW/REP1', # ok
			# 'K_HW/REP2', # ok
			# 'K_HW/REP3', # ok
			# 'K_HW/REP4', # ok
			# 'K_HW/REP5', # ok

			# 'Li/REP0', # ok
			# 'Li/REP1', # ok
			# 'Li/REP2', # ok
			# 'Li/REP3', # ok
			# 'Li/REP4', # ok
			# 'Li/REP5', # ok

			# 'Li_HW/REP0', # ok
			# 'Li_HW/REP1', # ok
			# 'Li_HW/REP2', # ok
			# 'Li_HW/REP3', # ok
			# 'Li_HW/REP4', # ok
			# 'Li_HW/REP5', # ok

			# 'Na/REP0_K', # ok 
			# 'Na/REP1_K', # ok
			# 'Na/REP2_K', # ok
			# 'Na/REP3_K', # ok
			# 'Na/REP4_K', # ok
			# 'Na/REP5_K', # ok
			# 'Na/REP0_Li', # ok
			# 'Na/REP1_Li', # ok
			# 'Na/REP2_Li', # ok
			# 'Na/REP3_Li', # ok
			# 'Na/REP4_Li', # ok
			# 'Na/REP5_Li', # ok

			# 'K/REP0_Na', # ok
			# 'K/REP1_Na', # not ok: 134 not converged
			# 'K/REP2_Na', # ok
			# 'K/REP3_Na', # ok
			# 'K/REP4_Na', # not ok: 141 not converged
			# 'K/REP5_Na', # ok
			# 'K/REP0_Li', # ok 
			# 'K/REP1_Li', # ok 
			# 'K/REP2_Li', # ok 
			# 'K/REP3_Li', # ok 
			# 'K/REP4_Li', # ok 
			# 'K/REP5_Li', # ok

			# 'Li/REP0_Na', # ok
			# 'Li/REP1_Na', # ok
			# 'Li/REP2_Na', # ok
			# 'Li/REP3_Na', # ok
			# 'Li/REP4_Na', # ok
			# 'Li/REP5_Na', # ok
			# 'Li/REP0_K', # ok 
			# 'Li/REP1_K', # ok 
			# 'Li/REP2_K', # ok 
			# 'Li/REP3_K', # not ok: forgot to submit it
			# 'Li/REP4_K', # ok 
			# 'Li/REP5_K', # ok 

			# 'Na_HW/REP0_K', # ok  
			# 'Na_HW/REP1_K', # ok 
			# 'Na_HW/REP2_K', # ok 
			# 'Na_HW/REP3_K', # ok 
			# 'Na_HW/REP4_K', # ok 
			# 'Na_HW/REP5_K', # ok 
			# 'Na_HW/REP0_Li', # ok 
			# 'Na_HW/REP1_Li', # ok 
			# 'Na_HW/REP2_Li', # ok 
			# 'Na_HW/REP3_Li', # ok 
			# 'Na_HW/REP4_Li', # ok
			# 'Na_HW/REP5_Li', # ok

			# 'K_HW/REP0_Na', # not ok: 131 not converged
			# 'K_HW/REP1_Na', # ok
			# 'K_HW/REP2_Na', # ok
			# 'K_HW/REP3_Na', # not ok: 108 not converged
			# 'K_HW/REP4_Na', # ok
			# 'K_HW/REP5_Na', # ok
			# 'K_HW/REP0_Li', # ok  
			# 'K_HW/REP1_Li', # ok  
			# 'K_HW/REP2_Li', # ok  
			# 'K_HW/REP3_Li', # ok  
			# 'K_HW/REP4_Li', # ok  
			# 'K_HW/REP5_Li', # ok 

			# 'Li_HW/REP0_K', # ok
			# 'Li_HW/REP1_K', # ok 
			# 'Li_HW/REP2_K', # ok 
			# 'Li_HW/REP3_K', # ok 
			# 'Li_HW/REP4_K', # ok 
			# 'Li_HW/REP5_K', # ok 
			# 'Li_HW/REP0_Na', # ok 
			# 'Li_HW/REP1_Na', # ok 
			# 'Li_HW/REP2_Na', # ok 
			# 'Li_HW/REP3_Na', # ok 
			# 'Li_HW/REP4_Na', # ok
			# 'Li_HW/REP5_Na', # ok

			# 'Na/REP0_B3LYP', # ok
			# 'Na_HW/REP0_B3LYP', # ok
			# 'K/REP0_B3LYP', # ok
			# 'K_HW/REP0_B3LYP', # ok
			# 'Li/REP0_B3LYP', # ok
			# 'Li_HW/REP0_B3LYP', # ok
			
			# 'Na/REP0_vaccum'
		]

		errors = []
		for key in keys:
			try:
				scan = scan_routine(key)
			except Exception as e:
				errors.append(f'{key}: {e}')
				mout.errorOut(f"Error encountered for {key}")
				mout.errorOut(e)
				continue
			pickle.dump(scan, open(f'{key}/scan.pickle','wb'))
			# break

		if errors:
			mout.headerOut("Errors:")
			for e in errors:
				mout.errorOut(e)

	else:
		# Select ref ion conformation
		ion = 'Na'

		pickle_files = sorted(glob.glob(f'{ion}/REP*/scan.pickle'))

		# remove B3LYP and vaccum files
		pickle_files = [f for f in pickle_files if len(f) < 23]

		Na_scans = []
		Na_K_scans = []
		Na_Li_scans = []

		for file in pickle_files:
			if file.split('/')[1][-1] == 'K':
				scan = pickle.load(open(file,'rb'))
				Na_K_scans.append(scan)
			elif file.split('/')[1][-2:] == 'Li':
				scan = pickle.load(open(file,'rb'))
				Na_Li_scans.append(scan)
			else:
				scan = pickle.load(open(file,'rb'))
				Na_scans.append(scan)

		atomic_rmsd(Na_scans[0], Na_scans[1])
		atomic_rmsd(Na_scans[0], Na_scans[1],align=False)
		exit()

		# Plot and store average heatmaps

		f1,f2,Na_avg = average_heatmaps(Na_scans,show=False)
		# mp.write(f'{ion}/average_heatmap.html',f1)
		# mp.write(f'{ion}/heatmap_rmsd.pdf',f2)

		f1,f2,Na_K_avg = average_heatmaps(Na_K_scans,show=False)
		# mp.write(f'{ion}/average_heatmap_k.html',f1)
		# mp.write(f'{ion}/heatmap_rmsd_k.pdf',f2)

		f1,f2,Na_Li_avg = average_heatmaps(Na_Li_scans,show=False)
		# mp.write(f'{ion}/average_heatmap_li.html',f1)
		# mp.write(f'{ion}/heatmap_rmsd_li.pdf',f2)

		subtract_heatmaps([Na_avg,Na_K_avg],show=True,ref_scan=Na_scans[0])

		exit()

		# RMSD between averages

		Na_K_rmsd = surface_rmsd(Na_avg, Na_K_avg)
		Na_Li_rmsd = surface_rmsd(Na_avg, Na_Li_avg)

		# RMSD ratios

		Na_ratio = Na_K_rmsd / Na_Li_rmsd
		mout.varOut(f"RMSD ratio between Na_K and Na_Li",Na_ratio)

		# Mean Difference between averages

		Na_K_md = surface_md(Na_avg, Na_K_avg)
		Na_Li_md = surface_md(Na_avg, Na_Li_avg)

		# Filtering average heatmaps
		min_threshold = -50
		max_threshold = 100
		print(f'Filtering data with min_threshold={min_threshold} kcal/mol and max_threshold={max_threshold} kcal/mol')

		Na_K_rmsd = surface_rmsd(Na_avg, Na_K_avg,threshold=True,min_threshold=min_threshold,max_threshold=max_threshold)
		Na_Li_rmsd = surface_rmsd(Na_avg, Na_Li_avg,threshold=True,min_threshold=min_threshold,max_threshold=max_threshold)
		Na_ratio = Na_K_rmsd / Na_Li_rmsd
		mout.varOut(f"RMSD ratio between Na_K and Na_Li",Na_ratio)

		Na_K_md = surface_md(Na_avg, Na_K_avg,threshold=True,min_threshold=min_threshold,max_threshold=max_threshold)
		Na_Li_md = surface_md(Na_avg, Na_Li_avg,threshold=True,min_threshold=min_threshold,max_threshold=max_threshold)
		Na_ratio = Na_K_md / Na_Li_md
		mout.varOut(f"MD ratio between Na_K and Na_Li",Na_ratio)

		exit()

		# Li_HW_data = compare_minima([
		# 	Li_HW_REP0,
		# 	Li_HW_REP1,
		# 	Li_HW_REP2,
		# 	Li_HW_REP3,
		# 	Li_HW_REP4,
		# 	Li_HW_REP5,
		# ])

		# compare_minima_energies([Na_data,Na_HW_data,K_data,K_HW_data,Li_data,Li_HW_data])

def compare_minima_energies(datasets,show=False):

	fig = go.Figure()

	names = [d['name'] for d in datasets]
	energy_diffs = [d['minima']['mean_z'] for d in datasets]
	energy_diffs_std = [d['minima']['std_z'] for d in datasets]

	trace = go.Bar(x=names,y=energy_diffs,error_y=dict(type='data',array=energy_diffs_std,visible=True))
	fig.add_trace(trace)

	fig.update_layout(yaxis_title='Relative minima depth kcal/mol')

	if show:
		fig.show()

	mp.write("minima_energies.pdf",fig)

def compare_minima(scans,show=False,calc_only=False):

	results = {}

	key = scans[0].key.split("/")[0]

	# ion and minima positions

	fig = go.Figure()

	minima_positions = []
	ion_positions = []

	for j,scan in enumerate(scans):

		x,y = scan.transform(scan.minimum.rcs)
		minima_positions.append([x,y])
		x,y = scan.transform(scan.reference_position)
		ion_positions.append([x,y])

		if j == 0:
			for i,(p,z) in enumerate(scan.pes_slices[:3]):

				p = scan.transform(p)

				x = [v[0] for v in p]
				y = [v[1] for v in p]

				trace = go.Scatter(name=f'{scan.key} slice {i}',legendgroup='slices',x=x,y=y,mode='lines',line=dict(color='black'))
				fig.add_trace(trace)

	ion_mean_x = np.mean([v[0] for v in ion_positions])
	ion_mean_y = np.mean([v[1] for v in ion_positions])
	ion_std_x = np.std([v[0] for v in ion_positions])
	ion_std_y = np.std([v[1] for v in ion_positions])

	min_mean_x = np.mean([v[0] for v in minima_positions])
	min_mean_y = np.mean([v[1] for v in minima_positions])
	min_std_x = np.std([v[0] for v in minima_positions])
	min_std_y = np.std([v[1] for v in minima_positions])

	if calc_only:
		return dict(min_mean_x=min_mean_x, min_mean_y=min_mean_y, min_std_x=min_std_x, min_std_y=min_std_y)

	fig.add_trace(go.Scatter(name='reference',mode='markers',x=[v[0] for v in ion_positions],y=[v[1] for v in ion_positions]))
	fig.add_trace(go.Scatter(name='minima',mode='markers',x=[v[0] for v in minima_positions],y=[v[1] for v in minima_positions]))
	
	fig.add_trace(go.Scatter(name='mean_reference',mode='markers',x=[ion_mean_x],y=[ion_mean_y],error_y=dict(type='data',array=[ion_std_x]),error_x=dict(type='data',array=[ion_std_y])))
	fig.add_trace(go.Scatter(name='mean_minima',mode='markers',x=[min_mean_x],y=[min_mean_y],error_y=dict(type='data',array=[min_std_x]),error_x=dict(type='data',array=[min_std_y])))

	fig.update_yaxes(scaleanchor="x",scaleratio=1)

	mp.write(f'{key}/minima_positions.pdf',fig)

	if show:
		fig.show()

	# relative minima depth
	energy_diffs = []
	minima_averages = []
	for scan in scans:
		minima = [scan.H0_H1_slice_min,scan.H0_H2_slice_min,scan.H1_H2_slice_min]
		avg_mins = np.mean(minima)
		minima_averages.append(avg_mins)
		delta = scan.interpolator(*scan.minimum.rcs) - avg_mins
		energy_diffs.append(delta)
	
	mean_energy_diff = np.mean(energy_diffs)
	std_energy_diff = np.std(energy_diffs)

	mean_minima_averages = np.mean(minima_averages)
	std_minima_averages = np.std(minima_averages)

	# area of surface within X kcal/mol of minima?

	# store results in the dictionary

	results['name'] = key
	results['ion'] = dict(mean_x=ion_mean_x,mean_y=ion_mean_y,std_x=ion_std_x,std_y=ion_std_y)
	results['minima'] = dict(mean_x=ion_mean_x,mean_y=ion_mean_y,std_x=ion_std_x,std_y=ion_std_y,mean_z=mean_energy_diff,std_z=std_energy_diff,mean_minima_averages=mean_minima_averages,std_minima_averages=std_minima_averages)

	return results

def average_heatmaps(scans,show=False,samples=50,plot_atoms=True):

	key = scans[0].key.split("/")[0]

	# get the transformed interpolator with absolute energies for each scan
	interpolators = [s.get_interpolator(transformed=True,absolute=True) for s in scans]

	# sampling space
	locut = min([s.transformed_locut for s in scans])
	hicut = max([s.transformed_hicut for s in scans])
	xv = np.linspace(locut, hicut, samples)
	yv = np.linspace(locut, hicut, samples)

	# sample points in the space and find the mean energies
	mean_points = []
	for y in yv:
		for x in xv:
			energies = [i(x,y) for i in interpolators]
			energies = [e for e in energies if not np.isnan(e)]
			if len(energies) == len(scans):
				mean_points.append([x,y,np.mean(energies)])

	locut = min([p[0] for p in mean_points])
	hicut = max([p[0] for p in mean_points])

	# generate a new interpolator from the mean points
	avg_interpolator = CloughTocher2DInterpolator([[p[0], p[1]] for p in mean_points],[p[2] for p in mean_points])

	# get the mean surface
	nx, ny = (100, 100)
	xv = np.linspace(min([p[0] for p in mean_points]), max([p[0] for p in mean_points]), nx)
	yv = np.linspace(min([p[1] for p in mean_points]), max([p[1] for p in mean_points]), ny)

	mesh_data = []
	for y in yv:

		y_slice = []
		for x in xv:
			energy = avg_interpolator(x,y)
			y_slice.append(energy)

		mesh_data.append(y_slice)

	mesh_data = np.array(mesh_data)

	fig = go.Figure()

	# average minima position
	data = compare_minima(scans,calc_only=True)
	fig.add_trace(go.Scatter(name='mean_minima',mode='markers',x=[data['min_mean_x']],y=[data['min_mean_y']],error_y=dict(type='data',array=[data['min_std_x']]),error_x=dict(type='data',array=[data['min_std_y']])))

	# plot the heatmap
	contours = dict(start=min([p[2] for p in mean_points]),end=max([p[2] for p in mean_points]),size=(max([p[2] for p in mean_points])-min([p[2] for p in mean_points]))/20)
	trace = go.Contour(x=xv,y=yv,z=mesh_data,contours=contours,showscale=False)
	fig.add_trace(trace)

	# plot the triangle
	trace = go.Scatter(x=[0,np.sqrt(3)/2,np.sqrt(3)/2,0],y=[0,0.5,-0.5,0],line=dict(color='black'),mode='lines')
	fig.add_trace(trace)

	fig.update_layout(title=f'average: {[s.key for s in scans]}')
	
	# add plot of atoms:
	if plot_atoms:
		fig = scans[0].example_sys.plot3d(fig=fig,flat=True,show=False,transform=scans[0].transform)

	if show:
		fig.show()

	# second figure plotting the RMSD between different surfaces

	fig2 = go.Figure()

	trace = go.Bar(x=[s.key for s in scans],y=[surface_rmsd(i, avg_interpolator,verbosity=0) for i in interpolators])
	fig2.add_trace(trace)

	fig2.update_layout(yaxis_title='SF PES RMSD [kcal/mol]')

	if show:
		fig2.show()

	return fig,fig2,avg_interpolator

def subtract_heatmaps(interpolators,show=False,samples=50,ref_scan=None,locut=-0.3,hicut=1.0):

	assert len(interpolators) == 2

	# key = scans[0].key.split("/")[0]

	# get the transformed interpolator with absolute energies for each scan
	# interpolators = [s.get_interpolator(transformed=True,absolute=True) for s in scans]

	# sampling space
	if ref_scan:
		locut = ref_scan.transformed_locut
		hicut = ref_scan.transformed_hicut
	# locut = min([s.transformed_locut for s in scans])
	# hicut = max([s.transformed_hicut for s in scans])
	xv = np.linspace(locut, hicut, samples)
	yv = np.linspace(locut, hicut, samples)

	# sample points in the space and find the subtracted energies
	sub_points = []
	for y in yv:
		for x in xv:
			energies = [i(x,y) for i in interpolators]
			energies = [e for e in energies if not np.isnan(e)]
			if len(energies) == 2:
				sub_points.append([x,y,energies[0]-energies[1]])

	locut = min([p[0] for p in sub_points])
	hicut = max([p[0] for p in sub_points])

	# generate a new interpolator from the subtracted points
	sub_interpolator = CloughTocher2DInterpolator([[p[0], p[1]] for p in sub_points],[p[2] for p in sub_points])

	# get the subtracted surface
	nx, ny = (100, 100)
	xv = np.linspace(min([p[0] for p in sub_points]), max([p[0] for p in sub_points]), nx)
	yv = np.linspace(min([p[1] for p in sub_points]), max([p[1] for p in sub_points]), ny)

	mesh_data = []
	for y in yv:

		y_slice = []
		for x in xv:
			energy = sub_interpolator(x,y)
			y_slice.append(energy)

		mesh_data.append(y_slice)

	mesh_data = np.array(mesh_data)

	fig = go.Figure()

	# # average minima position
	# data = compare_minima(scans,calc_only=True)
	# fig.add_trace(go.Scatter(name='mean_minima',mode='markers',x=[data['min_mean_x']],y=[data['min_mean_y']],error_y=dict(type='data',array=[data['min_std_x']]),error_x=dict(type='data',array=[data['min_std_y']])))

	# plot the heatmap
	contours = dict(start=min([p[2] for p in sub_points]),end=max([p[2] for p in sub_points]),size=(max([p[2] for p in sub_points])-min([p[2] for p in sub_points]))/20)
	trace = go.Contour(x=xv,y=yv,z=mesh_data,contours=contours,showscale=False)
	fig.add_trace(trace)

	# plot the triangle
	trace = go.Scatter(x=[0,np.sqrt(3)/2,np.sqrt(3)/2,0],y=[0,0.5,-0.5,0],line=dict(color='black'),mode='lines')
	fig.add_trace(trace)

	# fig.update_layout(title=f'average: {[s.key for s in scans]}')
	
	# add plot of atoms:
	if ref_scan:
		fig = ref_scan.example_sys.plot3d(fig=fig,flat=True,show=False,transform=ref_scan.transform)

	if show:
		fig.show()

	return fig,sub_interpolator

def surface_rmsd(interpolator1,interpolator2,x_range=[0,1],y_range=[-0.5,0.5],samples=50,verbosity=1,normalise=[np.sqrt(3)/4,0],threshold=False,min_threshold=-1,max_threshold=50):

	interpolators = [interpolator1,interpolator2]

	if normalise:
		references = [i(*normalise) for i in interpolators]

	# sampling space
	xv = np.linspace(x_range[0], x_range[1], samples)
	yv = np.linspace(y_range[0], y_range[1], samples)

	# sample points in the space and find the mean energies
	deltas = []
	for y in yv:
		for x in xv:

			energies = [i(x,y) for i in interpolators]

			energies = [e for e in energies if not np.isnan(e)]

			if normalise:
				energies = [e-r for e,r in zip(energies,references)]

			if threshold:
				energies = [e for e in energies if min_threshold < e < max_threshold]
							
			if len(energies) == 2:
				deltas.append(np.power(energies[0]-energies[1],2))

	rmsd = np.sqrt(np.mean(np.array(deltas)))

	if verbosity > 0:
		mout.varOut(f"RMSD from {len(deltas)} samples",rmsd,unit='kcal/mol')

	return rmsd

def surface_md(interpolator1,interpolator2,x_range=[0,1],y_range=[-0.5,0.5],samples=50,verbosity=1,normalise=[np.sqrt(3)/4,0],threshold=False,min_threshold=-1,max_threshold=50):

	interpolators = [interpolator1,interpolator2]

	if normalise:
		references = [i(*normalise) for i in interpolators]

	# sampling space
	xv = np.linspace(x_range[0], x_range[1], samples)
	yv = np.linspace(y_range[0], y_range[1], samples)

	# sample points in the space and find the mean energies
	deltas = []
	for y in yv:
		for x in xv:

			energies = [i(x,y) for i in interpolators]

			energies = [e for e in energies if not np.isnan(e)]

			if normalise:
				energies = [e-r for e,r in zip(energies,references)]
			
			if threshold:
				energies = [e for e in energies if min_threshold < e < max_threshold]
			
			if len(energies) == 2:
				deltas.append(energies[0]-energies[1])

	md = np.mean(np.array(deltas))

	if verbosity > 0:
		mout.varOut(f"Mean Difference from {len(deltas)} samples",md,unit='kcal/mol')

	return md

def atomic_rmsd(scan1,scan2,align=True):
	sys1 = scan1.example_sys.copy()
	sys2 = scan2.example_sys.copy()
	align_str = ''
	if align:
		align_str = '(aligned) '
		sys1.align_to(sys2)
	rmsd = sys1.rmsd(sys2)
	mout.varOut(f'Atomic RMSD {align_str}between {scan1.key} and {scan2.key}',rmsd,unit='Å')
	return rmsd

def scan_routine(key):

	scan = Scan(key)

	scan.get_reference_position()

	fig = scan.plot_heatmap(show_samples=True)
	mp.write(f'{scan.key}/heatmap.html',fig)

	scan.minimum = hijack.FakeAtoms(rcs=[scan.reference_position[0],scan.reference_position[1]],name="minimum")

	scan.optimize(scan.minimum)

	interpolator = scan.get_interpolator()

	mout.varOut(f"Optimised ({scan.minimum.name}) RCS",scan.minimum.rcs)
	mout.varOut(f"Optimised ({scan.minimum.name}) energy",interpolator(*scan.minimum.rcs))

	scan.plot_energy_slices(show=False)
	
	fig = scan.plot_heatmap(show=False,show_missing=True)

	mp.write(f'{scan.key}/heatmap_min.html',fig)

	fig = scan.plot_heatmap(show=False,transformed=True,show_samples=True)
	mp.write(f'{scan.key}/heatmap_transformed.html',fig)

	# exit()

	return scan

if __name__ == '__main__':
	main()

