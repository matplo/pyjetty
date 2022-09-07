from __future__ import absolute_import, division, print_function
import sys
import os
import argparse
import fastjet as fj
import math
import numpy as np
import math
from matplotlib.path import Path
from matplotlib.patches import PathPatch
# from scipy.spatial import ConvexHull
# from scipy.spatial import KDTree
import matplotlib.lines as mlines
import matplotlib.pyplot as plt


def circle_line_segment_intersection(circle_center, circle_radius, pt1, pt2, full_line=True, tangent_tol=1e-9):
    """ Find the points at which a circle intersects a line-segment.  This can happen at 0, 1, or 2 points.

    :param circle_center: The (x, y) location of the circle center
    :param circle_radius: The radius of the circle
    :param pt1: The (x, y) location of the first point of the segment
    :param pt2: The (x, y) location of the second point of the segment
    :param full_line: True to find intersections along full line - not just in the segment.  False will just return intersections within the segment.
    :param tangent_tol: Numerical tolerance at which we decide the intersections are close enough to consider it a tangent
    :return Sequence[Tuple[float, float]]: A list of length 0, 1, or 2, where each element is a point at which the circle intercepts a line segment.

    Note: We follow: http://mathworld.wolfram.com/Circle-LineIntersection.html
    """

    (p1x, p1y), (p2x, p2y), (cx, cy) = pt1, pt2, circle_center
    (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
    dx, dy = (x2 - x1), (y2 - y1)
    dr = (dx ** 2 + dy ** 2)**.5
    big_d = x1 * y2 - x2 * y1
    discriminant = circle_radius ** 2 * dr ** 2 - big_d ** 2

    if discriminant < 0:  # No intersection between circle and line
        return []
    else:  # There may be 0, 1, or 2 intersections with the segment
        intersections = [
            (cx + (big_d * dy + sign * (-1 if dy < 0 else 1) * dx * discriminant**.5) / dr ** 2,
             cy + (-big_d * dx + sign * abs(dy) * discriminant**.5) / dr ** 2)
            for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
        if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
            fraction_along_segment = [
                (xi - p1x) / dx if abs(dx) > abs(dy) else (yi - p1y) / dy for xi, yi in intersections]
            intersections = [pt for pt, frac in zip(
                intersections, fraction_along_segment) if 0 <= frac <= 1]
        # If line is tangent to circle, return just one point (as both intersections have same location)
        if len(intersections) == 2 and abs(discriminant) <= tangent_tol:
            return [intersections[0]]
        else:
            return intersections


class SubjetPatchPlot(object):
	def __init__(self, xy, radius, clipR, npoints=20, name='SubjetPatchPlot'):
		self.name = name
		self.radius = radius
		self.xy = (xy[0], xy[1])
		self.center = self.xy
		self.npoints = npoints
		self.clipR = clipR
		self.clipped_to = []
		d = math.dist(self.xy, [0., 0.])
		if d > clipR:
			print('[error] center outside the clipping area?')
			return False

	def clip_to_R(self, x, y, clipR):
		d = math.dist((x, y), (0.0, 0.0))
		if d > clipR:
			isec = circle_line_segment_intersection((0, 0), clipR, self.center, (x, y))
			if len(isec) < 2:
				print('strange: point not intersecting twice with jet R')
				return x, y
			_d1 = math.dist((isec[0][0], isec[0][1]), self.center)
			_d2 = math.dist((isec[1][0], isec[1][1]), self.center)
			if _d1 < _d2:
				x, y = isec[0][0], isec[0][1]
			else:
				x, y = isec[1][0], isec[1][1]
		return x, y

	def make_circle_points(self, center, r, npoints=-1):
		if npoints < 0:
			npoints = self.npoints
		t = np.linspace(0, 2. * np.pi, npoints)
		_tmp_x = []
		_tmp_y = []
		_tmp_t = []
		for _t in t:
			_x = center[0] + r * math.cos(_t)
			_y = center[1] + r * math.sin(_t)
			# _x, _y = self.clip_to_R(_x, _y, self.clipR)
			if math.dist((_x, _y), (0.0, 0.0)) > self.clipR:
				_x = self.clipR * math.cos(_t)
				_y = self.clipR * math.sin(_t)
				if math.dist((_x, _y), center) > r:
					continue
			_tmp_x.append(_x)
			_tmp_y.append(_y)
			_tmp_t.append(_t)
		return _tmp_x, _tmp_y, _tmp_t

	def make_points(self):
		self.t_orig = np.linspace(0, 2. * np.pi, self.npoints)
		self.x = []
		self.y = []
		self.t = []

		self.x, self.y, self.t = self.make_circle_points(
			self.center, self.radius)
		self.vertices = np.array([self.x, self.y]).T
		# self.hull = ConvexHull(self.vertices)
		# self.kdtree = KDTree(self.vertices)
		self.codes = np.full(len(self.vertices), Path.LINETO)
		self.codes[len(self.vertices)-1] = Path.MOVETO
		self.codes[0] = Path.MOVETO

	def plot_patch_ax(self, ax, color=[0.3, 0, 0.0, 0.5], points=False, line_width=0):
		self.make_points()
		self.path = Path(self.vertices, self.codes)
		self.path_patch = PathPatch(self.path, facecolor=color, lw=line_width)
		ax.add_patch(self.path_patch)
		if points:
			ax.plot(self.x, self.y, 'p')
			# ax.plot(self.vertices[self.hull.vertices, 0], self.vertices[self.hull.vertices, 1], 'r--', lw=2)


class SingleSubjetPlot(object):
	def __init__(self, j, sj, sj_r, color, clipR=0.4, npoints=50, name='SingleSubjetPlot'):
		self.name = name
		self.npoints = npoints
		self.color = color
		self.clipR = clipR
		self.sj = sj
		self.sj_r = sj_r
		self.j = j
		self.sj_dphi = self.j.delta_phi_to(self.sj)
		self.sj_deta = self.j.rapidity() - self.sj.rapidity()
		self.sjpatch = SubjetPatchPlot([self.sj_dphi, self.sj_deta], self.sj_r, self.clipR, self.npoints, name=self.name+'Patch')

	def plot_subjet(self, ax, scale_pt = 0):
		self.pts = []
		self.ys = []
		self.phis = []
		self.cs = []
		self.lines = []
		self.circles = []
		self.sc = fj.sorted_by_pt(self.sj.constituents())
		self.pts.extend([p.perp() for p in self.sc])
		self.ys.extend([self.j.rapidity() - p.rapidity() for p in self.sc])
		self.phis.extend([self.j.delta_phi_to(p) for p in self.sc])
		for p in self.sc:
			self.cs.append(self.color)
			part_dphi = self.j.delta_phi_to(p)
			part_deta = self.j.rapidity() - p.rapidity()
			self.lines.append([self.sj_dphi, part_dphi,
							self.sj_deta, part_deta,
							self.cs[-1]])

		self.circles.append([self.sj_dphi, self.sj_deta, self.sj_r, self.color])

		if scale_pt == 0:
			scale_pt = self.j.perp()
		self.zs = [pt/scale_pt for pt in self.pts]
		self.zs_sized = [z*1000. for z in self.zs]
		ax.scatter(self.phis, self.ys, c=self.cs, s=self.zs_sized, alpha=0.4, cmap="magma")  # cmap='viridis')

		transform = ax.transData
		for l in self.lines:
			_line = mlines.Line2D([l[0], l[1]], [l[2], l[3]], color=l[4])
			_line.set_transform(transform)
			ax.add_line(_line)
	
	def plot_patch(self, ax, color=None):
		if color is None:
			color = self.color
		self.sjpatch.plot_patch_ax(ax, color=color, points=False, line_width=0)



class SubjetPlot(object):
	def __init__(self, jet, clipR, npoints, bg_color=[1.0, 1.0, 1.0, 1.0]):
		self.jet = jet
		self.subjets = []
		self.clipR = clipR
		self.npoints = npoints
		self.bg_color = bg_color
  
	def set_colors(self):
		self.colors = [[1, 0, 0, 0.3], [0.1, .75, 0.1, 0.3], [0, 0, 1, 0.3]]
		for i in range(3, len(self.subjets)):
			gr_col = 0.1 * (i - 2)
			if gr_col > 0.9:
				gr_col = 0.9
			_col = [gr_col, gr_col, gr_col, 0.3]
			print(i, _col)
			self.colors.append(_col)

	def plot(self, scale_pt=0, sj_r=0.1):

		self.subjets = []
		sj_def = fj.JetDefinition(fj.antikt_algorithm, sj_r)
		self.subjets = fj.sorted_by_pt(sj_def(self.jet.constituents()))
		self.set_colors()

		jet_R0 = self.clipR
		phis = []
		phis.append(jet_R0)
		phis.append(-jet_R0)
		ys = []
		ys.append(jet_R0)
		ys.append(-jet_R0)
		pts = []
		pts.append(0)
		pts.append(0)
		if scale_pt == 0:
			scale_pt = self.jet.perp()
		zs = [pt/scale_pt for pt in pts]
		zs_sized = [z*1000. for z in zs]
		cs = []
		cs.append([0,0,0,0])
		cs.append([0,0,0,0])

		self.fig = plt.figure()
		plt.scatter(phis, ys, c=cs, s=zs_sized, alpha=0.4, cmap="magma") #cmap='viridis')

		plt.xlabel(r"$\Delta\varphi$")
		plt.ylabel(r"$\Delta y$")

		ax = self.fig.axes[0]

		self.sj_subplots = []
		for i,sj in enumerate(self.subjets):
			_sjplt = SingleSubjetPlot(self.jet, sj, sj_r, self.colors[i], self.clipR, self.npoints, name='sjplt{}'.format(i))
			self.sj_subplots.append(_sjplt)

		for i,sj in enumerate(reversed(self.sj_subplots)):
			sj.plot_patch(ax)
			for j, sjB in enumerate(reversed(self.sj_subplots)):
				if j <= i:
					continue
				sjB.plot_patch(ax, color=self.bg_color)
			sj.plot_subjet(ax)

		ax.set_xlim(-0.5, 0.5)
		ax.set_ylim(-0.5, 0.5)
		ax.set(aspect=1)
