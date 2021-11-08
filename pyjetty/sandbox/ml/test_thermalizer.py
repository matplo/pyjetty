#!/usr/bin/env python

import fastjet as fj

import ROOT
ROOT.gSystem.Load('libpyjetty_rutil')

# standard numerical library imports
import numpy as np
rng = np.random.RandomState(0)

# matplotlib is required for this example
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10,9)
from matplotlib import animation, rc

from IPython.core.display import display, HTML
display(HTML("<style>div.output_scroll { height: 200em; }</style>"))


class scatter():
	def __init__(self,x,y,ax,size=1,**kwargs):
		self.n = len(x)
		self.ax = ax
		self.ax.figure.canvas.draw()
		self.size_data=size
		self.size = size
		self.sc = ax.scatter(x,y,s=self.size,**kwargs)
		self._resize()
		self.cid = ax.figure.canvas.mpl_connect('draw_event', self._resize)

	def _resize(self,event=None):
		ppd=72./self.ax.figure.dpi
		trans = self.ax.transData.transform
		s =  ((trans((1,self.size_data))-trans((0,0)))*ppd)[1]
		if s != self.size:
			self.sc.set_sizes(s**2*np.ones(self.n))
			self.size = s
			self._redraw_later()
	
	def _redraw_later(self):
		self.timer = self.ax.figure.canvas.new_timer(interval=10)
		self.timer.single_shot = True
		self.timer.add_callback(lambda : self.ax.figure.canvas.draw_idle())
		#self.timer.start()


def draw_jet(j, Rabs=0.4):
	# for the jet but not subjets
	pts = [p.perp() for p in fj.sorted_by_pt(j.constituents())]
	ys = [j.rapidity() - p.rapidity() for p in fj.sorted_by_pt(j.constituents())]
	phis = [j.delta_phi_to(p) for p in fj.sorted_by_pt(j.constituents())]

	print("n | pt | z | phi | eta | dR")
	for i, p in enumerate(fj.sorted_by_pt(j.constituents())):
		print(i, p.perp(), p.perp()/j.perp(), p.phi(), p.eta(), p.delta_R(j))
	
	phis.append(Rabs)
	phis.append(-Rabs)
	ys.append(Rabs)
	ys.append(-Rabs)
	pts.append(0)
	pts.append(0)

	zs 			= [pt/j.perp() for pt in pts]
	zs_sized 	= [z * 1000. for z in zs]
	cs 			= [int(z * 100.) for z in zs]

	fig = plt.figure()

	# plt.scatter(phis, ys, c=colors, s=zs_sized, alpha=0.4, cmap="PuOr") #cmap='viridis')
	# plt.scatter(phis, ys, c=cs, s=zs_sized, alpha=0.4, cmap="magma") #cmap='viridis')
	plt.scatter(phis, ys, c=cs, s=zs_sized, alpha=0.4, cmap='viridis')
	plt.xlabel('phi')
	plt.ylabel('y')
	# plt.colorbar();  # show color scale
	plt.show()


def draw_parts_thermal(parts, Rabs=0.4):
	jet_R0 = Rabs * 1.0
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	parts_v = [fj.PseudoJet(p.px(), p.py(), p.pz(), p.e()) for p in parts]
	jets = jet_def(parts_v)
	for j in jets:
		draw_jet(j, jet_R0)
		break


th0 = ROOT.RUtil.Thermalizer()
parts = th0.thermalize(10, 0, 0)
draw_parts_thermal(parts, 10)

th1 = ROOT.RUtil.Thermalizer(0.7, -1, 0.4, -1)
parts = th1.thermalize(10, 0, 0)
draw_parts_thermal(parts, 0.4)

thM = ROOT.RUtil.Thermalizer(0.7, 3, 0.4, -1)
parts = thM.thermalize(10, 0, 0)
draw_parts_thermal(parts, 0.4)

thM = ROOT.RUtil.Thermalizer(0.7, 5, 0.4, -1)
parts = thM.thermalize(10, 0, 0)
draw_parts_thermal(parts, 0.4)

thM = ROOT.RUtil.Thermalizer(0.7, 2, 0.4, -1)
parts = thM.thermalize(10, 0, 0)
draw_parts_thermal(parts, 0.4)
