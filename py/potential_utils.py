# encoding: utf-8
# 2022 © Vasileios Angelidakis <vasileios.angelidakis@ncl.ac.uk>
"""
Auxiliary functions for the Potential Blocks
"""

import math, random, doctest, geom, numpy
from yade import Vector3, Quaternion, utils
from yade.wrapper import *
from math import sin, cos, tan, sqrt, pi, radians  # atan, atan2

import logging

logging.basicConfig(level=logging.INFO)

from numpy import array
#from yade import utils	# FIXME: I only use utils._commonBodySetup
from yade.utils import _commonBodySetup, randomColor

try:  # use psyco if available
	import psyco
	psyco.full()
except ImportError:
	pass


#**********************************************************************************
#creates Potential Blocks, defining the coefficients of their faces
def potentialblock(material, a=[], b=[], c=[], d=[], r=0.0, R=0.0, mask=1, isBoundary=False, fixed=False, color=[-1, -1, -1]):
	"""creates potential block.

	:param Material material: material of new body
	:param [float] a,b,c,d: lists of plane coefficients of the particle faces (see PotentialBlock docs)
	:param float r: radius of inner Potential Particle (see PotentialBlock docs)
	:param float R: distance R of the Potential Blocks (see PotentialBlock docs)
	:param bool isBoundary: whether this is a boundary body (see PotentialBlock docs)
	"""
	# TODO: In this function, I can introduce all the other attributes of the PBs, like: fixedNormal, boundaryNormal -> Naah, the invocation gets too long. I can use *kw!!!
	pb = Body()
	pb.mask = mask
	pb.aspherical = True

	pb.shape = PotentialBlock(
	        a=a, b=b, c=c, d=d, r=r, R=R, isBoundary=isBoundary, AabbMinMax=True
	)  #id=len(O.bodies) #FIXME: Check if I need id for vtk output
	if color[0] == -1:
		pb.shape.color = randomColor(seed=random.randint(0, 1E6))
	else:
		pb.shape.color = color
	utils._commonBodySetup(pb, pb.shape.volume, pb.shape.inertia, material, pos=pb.shape.position, fixed=fixed)
	pb.state.ori = pb.shape.orientation
	return pb


#**********************************************************************************
#creates cuboidal particle using the Potential Blocks
def cuboid(material, edges=Vector3(1, 1, 1), r=0.0, R=0.0, center=[0, 0, 0], mask=1, isBoundary=False, fixed=False, color=[-1, -1, -1]):
	"""creates cuboid using the Potential Blocks

	:param Vector3 edges: edges of the cuboid
	:param Material material: material of new body (FrictMat)
	:param Vector3 center: center of the new body
	"""
	aa = [1, -1, 0, 0, 0, 0]
	bb = [0, 0, 1, -1, 0, 0]
	cc = [0, 0, 0, 0, 1, -1]
	dd = [edges[0] / 2., edges[0] / 2., edges[1] / 2., edges[1] / 2., edges[2] / 2., edges[2] / 2.]

	if not r:
		r = min(dd) / 2.
	cuboid = potentialblock(material=material, a=aa, b=bb, c=cc, d=array(dd) - r, r=r, R=R, mask=mask, isBoundary=isBoundary, fixed=fixed, color=color)
	cuboid.state.pos = center
	return cuboid


#**********************************************************************************
#creates Aabb boundary plates using the Potential Blocks, around a given cuboidal space
def aabbPlates(material, extrema=None, thickness=0.0, r=0.0, R=0.0, mask=1, isBoundary=False, fixed=True, color=None):
	"""Return 6 cuboids that will wrap existing packing as walls from all sides.			#FIXME: Correct this comment

	:param Material material: material of new bodies (FrictMat)
	:param [Vector3, Vector3] extrema: extremal points of the Aabb of the packing, as a list of two Vector3, or any equivalent type (will not be calculated if not specified)
	:param float thickness: wall thickness (equal to 1/10 of the smallest dimension if not specified)
	:param float r: radius of inner Potential Particle (see PotentialBlock docs)
	:param float R: distance R of the Potential Blocks (see PotentialBlock docs)
	:param int mask: groupMask for the new bodies

	:returns: a list of 6 PotentialBlock Bodies enclosing the packing, in the order minX,maxX,minY,maxY,minZ,maxZ.
	"""
	walls = []
	#	if not extrema: extrema=aabbExtrema() #TODO: aabbExtrema() is not compatible with non-spherical particles yet
	if not thickness:
		thickness = min(extrema[1][0] - extrema[0][0], extrema[1][1] - extrema[0][1], extrema[1][2] - extrema[0][2]) / 10.

	randColor = False
	if not color:
		randColor = True

	for axis in [0, 1, 2]:
		mi, ma = extrema
		center = [(mi[i] + ma[i]) / 2. for i in range(3)]
		extents = [(ma[i] - mi[i]) for i in range(3)]
		extents[axis] = thickness / 2.
		if randColor:
			color = randomColor(seed=random.randint(0, 1E6))
		for j in [0, 1]:
			center[axis] = extrema[j][axis] + (j - .5) * thickness / 2.
			walls.append(
			        cuboid(material=material, edges=extents, r=r, R=R, center=center, mask=mask, isBoundary=isBoundary, fixed=fixed, color=color)
			)
			walls[-1].shape.wire = True
	return walls


#**********************************************************************************
#creates cylindrical boundary plates using the Potential Blocks, using a given radius and around a given axis
def cylindricalPlates(
        material, radius=0.0, height=0.0, thickness=0.0, numFaces=3, r=0.0, R=0.0, mask=1, isBoundary=False, fixed=True, lid=[True, True], color=None
):
	"""Return numFaces cuboids that will wrap existing packing as walls from all sides. 			#FIXME: Correct this comment

	:param Material material: material of new bodies (FrictMat)
	:param float radius: radius of the cylinder
	:param float height: height of cylinder
	:param float thickness: thickness of cylinder faces (equal to 1/10 of the cylinder inradius if not specified)
	:param int numFaces: number of cylinder faces (>3)
	:param float r: radius of inner Potential Particle (see PotentialBlock docs)
	:param float R: distance R of the Potential Blocks (see PotentialBlock docs)
	:param int mask: groupMask for the new bodies
	:param lid [bool]: list of booleans, whether to create the bottom and top  lids of the cylindrical plates

	:returns: a list of cuboidal Potential Blocks forming a hollow cylinder
	"""
	# TODO: Check facetCylinder for orientation of the cylinder
	# TODO: Add center of the cylinder

	walls = []
	if not thickness:
		thickness = min(radius, height / 2.) / 10.
	angle = radians(360) / numFaces
	axis = Vector3(0, 0, 1)  #TODO: To make it work for any axis - have to change: center, edges

	randColor = False
	if not color:
		randColor = True
	for i in range(0, numFaces):
		center = Vector3((radius + thickness / 2.) * cos(i * angle), (radius + thickness / 2.) * sin(i * angle), height / 2.)  #*(axis.asDiagonal())

		if randColor:
			color = [abs(cos(i * angle)), abs(sin(i * angle)), 1]
		walls.append(
		        cuboid(
		                material=material,
		                edges=Vector3(thickness, 2 * (radius) * tan(pi / numFaces), height),
		                r=r,
		                R=R,
		                center=center,
		                mask=mask,
		                isBoundary=isBoundary,
		                fixed=fixed,
		                color=color
		        )
		)
		walls[-1].state.ori = Quaternion(axis, i * angle)
		walls[-1].shape.wire = True

	color = [0.36, 0.54, 0.66]
	if lid[0] == True:
		# Create top plate
		walls.append(
		        prism(
		                material=material,
		                radius1=radius,
		                thickness=thickness,
		                numFaces=numFaces,
		                r=r,
		                R=R,
		                color=color,
		                mask=mask,
		                isBoundary=isBoundary,
		                fixed=fixed
		        )
		)
		walls[-1].state.pos = axis * (height + thickness / 2)
		walls[-1].state.ori = Quaternion(axis, i * angle)  #FIXME: Here I use i outside the loop?
	if lid[1] == True:
		# Create bottom plate
		walls.append(
		        prism(
		                material=material,
		                radius1=radius,
		                thickness=thickness,
		                numFaces=numFaces,
		                r=r,
		                R=R,
		                color=color,
		                mask=mask,
		                isBoundary=isBoundary,
		                fixed=fixed
		        )
		)
		walls[-1].state.pos = -axis * thickness / 2
		walls[-1].state.ori = Quaternion(axis, i * angle)  #FIXME: Here I use i outside the loop?

	return walls


##**********************************************************************************
#creates regular prism with N faces
def prism(material, radius1=0.0, radius2=-1, thickness=0.0, numFaces=3, r=0.0, R=0.0, center=None, color=[1, 0, 0], mask=1, isBoundary=False, fixed=False):
	"""Return regular prism with numFaces

	:param Material material: material of new bodies (FrictMat)
	:param float radius1: inradius of the start cross-section of the prism
	:param float radius2: inradius of the end cross-section of the prism (equal to radius1 if not specified)
	:param float thickness: thickness of the prism (equal to radius1 if not specified)
	:param int numFaces: number of prisms' faces (>3)
	:param float r: radius of inner Potential Particle (see PotentialBlock docs)
	:param float R: distance R of the Potential Blocks (see PotentialBlock docs)
	:param Vector3 center: center of the Potential Blocks (not currently used)
	:param int mask: groupMask for the new bodies

	:returns: an axial-symmetric Potential Block with variable cross-section, which can become either a regular prism (radius1=radius2), a pyramid (radius2=0) or a cylinder or cone respectively, for a large enough numFaces value.
	"""
	aa = []
	bb = []
	cc = []
	dd = []
	if radius2 == -1:
		radius2 = radius1
	if not thickness:
		thickness = radius1
	angle = radians(360) / numFaces

	for i in range(0, numFaces):
		aTemp = cos(i * angle)
		bTemp = sin(i * angle)
		cTemp = (radius1 - radius2) / thickness
		dTemp = (radius1 + radius2) / 2

		magnitude = Vector3(aTemp, bTemp, cTemp).norm()

		aa.append(aTemp / magnitude)
		bb.append(bTemp / magnitude)
		cc.append(cTemp / magnitude)
		dd.append(dTemp / magnitude)

	aa.extend([0.0, 0.0])
	bb.extend([0.0, 0.0])
	cc.extend([1.0, -1.0])
	dd.extend([thickness / 2., thickness / 2.])

	if not r:
		r = min(dd) / 2.
	prism = potentialblock(material=material, a=aa, b=bb, c=cc, d=array(dd) - r, r=r, R=R, mask=mask, isBoundary=isBoundary, fixed=fixed, color=color)
	#	if center: prism.state.pos=prism.state.pos+center; print(center) #FIXME: Maybe assign center if not (0,0,0), but add it to the local position, rather than overwriting it
	return prism


##**********************************************************************************
#TODO: export PotentialBlock to stl file: This would better fit in the export module

##**********************************************************************************
#TODO: creates masonry arch: I can do it in two ways:
#	1. Geometrically, using the angles or
#	2. Using BlockGen, demonstrating that we can do non-permanent joints
#	Explore multi-ring arches: For the inner ring, I must replicate the interlocking geometry

##**********************************************************************************
#TODO: creates a model for periodic and aperiodic masonry

##**********************************************************************************
#TODO: creates SAG mill with chosen: radius, width, optionally rigid boundaries or no boundaries, dent size

##**********************************************************************************
#creates PotentialBlock or PotentialParticle from polyhedra # FIXME: Unfinished
#def polyhedra2PotentialBlock(b,r=None,R=0.0): #maybe polyhedra2Potential or polyhedra2PotentialShape #FIXME: Haven't checked it yet
#	"""Calculate coefficients of the planes comprising the faces of a polyhedron, in the format used by the PotentialBlock and PotentialParticle shape classes

#	:param Body b:			Body with b.shape=Polyhedra
#	:param vector(Real) aa:	Coefficients 'a' of planes: ax+by+cz-d-r=0
#	:param vector(Real) bb:	Coefficients 'b' of planes: ax+by+cz-d-r=0
#	:param vector(Real) cc:	Coefficients 'c' of planes: ax+by+cz-d-r=0
#	:param vector(Real) dd:	Coefficients 'd' of planes: ax+by+cz-d-r=0
#	:param Real r:			Suggested radius of PotentialBlock or PotentialParticle to be generated

#	"""
#	Faces=b.shape.GetSurfaces()
#	Vertices = [b.state.ori*v for v in b.shape.v] # vertices in global coords

#	aa=[]
#	bb=[]
#	cc=[]
#	dd=[]
#	for fv in range(len(Faces)):
#		v1=Vertices[Faces[fv][0]] #1st vertex on plane fv
#		v2=Vertices[Faces[fv][1]] #2nd vertex on plane fv
#		v3=Vertices[Faces[fv][2]] #3rd vertex on plane fv

#		n1=v2-v1 #1st vector on face fv
#		n2=v3-v1 #2nd vector on face fv
#		nNormal=n1.cross(n2) #normal vector of face fv
#		nNormal.normalize()

#		aTemp=nNormal[0]
#		bTemp=nNormal[1]
#		cTemp=nNormal[2]

#		dTemp = (aTemp*v1[0] + bTemp*v1[1] + cTemp*v1[2])

#		aa.append(aTemp)
#		bb.append(bTemp)
#		cc.append(cTemp)
#		dd.append(dTemp)

##	if not r: r=0.5*abs(min(dd)) # Recommended value. This assignment here is only meaningful if the particle is centered to its centroid

#	pb=potentialblock(material=material,a=aa,b=bb,c=cc,d=array(dd)-r,r=r,R=R,mask=mask,isBoundary=isBoundary,fixed=fixed,color=color)
#	return pb

##	b2=Body()
##	b2.aspherical=True
##	color=Vector3(random.random(),random.random(),random.random())
##	b2.shape=PotentialBlock(k=0.0, r=r, R=0.0, a=aa, b=bb, c=cc, d=array(dd)-r, AabbMinMax=True, color=color)
##	utils._commonBodySetup(b2, b2.shape.volume, b2.shape.inertia, material='frictmat', pos=b1.state.pos, fixed=False)
##	b2.state.ori=b2.shape.orientation
##	O.bodies.append(b2)
#	return aa,bb,cc,dd,r

##**********************************************************************************
#TODO: creates polyhedra from PotentialBlock or PotentialParticle

#from numpy import arange
##**********************************************************************************
##TODO: creates ellipsoidal particle using the Potential Blocks # FIXME: Unfinished
#def ellipsoid(material, axes=Vector3(1,1,1), numFaces=20, r=0.0, R=None, center=[0,0,0], mask=1, isBoundary=False, fixed=False, color=[-1,-1,-1]):
#	"""creates polyhedral ellipsoid using the Potential Blocks

#	:param Vector3 axes: axes of the ellipsoid
#	:param int numFaces: number of particle faces (the total number of faces will be near this value, not exactly it)
#	:param Material material: material of new body (FrictMat)
#	:param Vector3 center: center of the new body
#	"""
#	thetaStep=pi/8.
#	phiStep=pi/8.

#	aa=[]; bb=[]; cc=[]; dd=[];
##	for theta in arange(0, 2*pi, thetaStep):
##		for phi in arange(-pi/2+pi/4., pi/2., phiStep):
##			aa.extend([cos(phi)*s(theta)])
##			bb.extend([cos(phi)*sin(theta)])
##			cc.extend([sin(phi)])
##			dd.extend([1.-r])

#	for theta in arange(0, pi+thetaStep, thetaStep):
#		for phi in arange(0, 2*pi, phiStep):
#			aa.extend([sin(theta)*cos(phi)])
#			bb.extend([sin(theta)*sin(phi)])
#			cc.extend([cos(theta)])
#			dd.extend([1.-r])

#	# FIXME: Here, I have to correct the edge size, using a scaleFactor
#	# FIXME: I also have to fix the shape itself

#	ellipsoid = potentialblock(material=material,a=aa,b=bb,c=cc,d=dd,r=r,R=R,mask=mask,isBoundary=isBoundary,fixed=fixed,color=color)
#	ellipsoid.state.pos = center
#	return ellipsoid


#**********************************************************************************
#creates platonic solids using the Potential Blocks
def platonic_solid(
        material,
        numFaces,
        edge=0.0,
        ri=0.0,
        rm=0.0,
        rc=0.0,
        volume=0.0,
        r=0.0,
        R=None,
        center=[0, 0, 0],
        mask=1,
        isBoundary=False,
        fixed=False,
        color=[-1, -1, -1]
):
	errors = 0
	"""creates platonic solids using the Potential Blocks
	User must specify either the edge, the inradius, the circumradius or the volume of the particle (only one of them)

	:param int numFaces: number of particle faces (regular 4: tetrahedron, 6: hexahedron (cube), 8: octahedron, 12: dodecahedron, 20: icosahedron)
	:param float edge: edge of the platonic solid
	:param float ri: inradius of the platonic solid
	:param float rm: midradius of the platonic solid
	:param float rc: circumradius of the platonic solid
	:param float volume: volume of the platonic solid
	:param Material material: material of new body (FrictMat)
	:param Vector3 center: center of the new body
	"""
	inputParams = [edge, ri, rm, rc, volume]
	count = sum(1 for i in inputParams if i > 0.0)
	if (count > 1):
		logging.error(' Assign only one of: edge - ri - rm - rc - volume')
		return (None)

	gamma = 1 / sqrt(3)

	delta = sqrt((5 - sqrt(5)) / 10.)
	epsilon = sqrt((5 + sqrt(5)) / 10.)

	zeta = sqrt((3 - sqrt(5)) / 6.)
	eta = sqrt((3 + sqrt(5)) / 6.)

	# Schläfli symbols {p,q}: https://en.wikipedia.org/wiki/Platonic_solid

	p = {4: 3, 6: 4, 8: 3, 12: 5, 20: 3}
	q = {4: 3, 6: 3, 8: 4, 12: 3, 20: 5}
	h = {4: 4, 6: 6, 8: 6, 12: 10, 20: 10}  # Coxeter number
	t = cos(pi / q[numFaces]) / sin(pi / h[numFaces])  # tan(theta/2)

	if (numFaces == 4):  #tetrahedron
		aa = array([+gamma, +gamma, -gamma, -gamma])
		bb = array([-gamma, +gamma, +gamma, -gamma])
		cc = array([+gamma, -gamma, +gamma, -gamma])

	elif (numFaces == 6):  #hexahedron (cube)
		aa = [1, -1, 0, 0, 0, 0]
		bb = [0, 0, 1, -1, 0, 0]
		cc = [0, 0, 0, 0, 1, -1]

	elif (numFaces == 8):  #octahedron
		aa = array([+gamma, -gamma, -gamma, +gamma, +gamma, -gamma, +gamma, -gamma])
		bb = array([+gamma, -gamma, +gamma, -gamma, -gamma, +gamma, +gamma, -gamma])
		cc = array([+gamma, -gamma, +gamma, -gamma, +gamma, -gamma, -gamma, +gamma])

	elif (numFaces == 12):  #dodecahedron
		aa = array([delta, -delta, delta, -delta, 0, 0, 0, 0, epsilon, -epsilon, -epsilon, epsilon])
		bb = array([epsilon, -epsilon, -epsilon, epsilon, delta, -delta, delta, -delta, 0, 0, 0, 0])
		cc = array([0, 0, 0, 0, epsilon, -epsilon, -epsilon, epsilon, delta, -delta, delta, -delta])

	elif (numFaces == 20):  #icosahedron: first 8 faces
		aa = [
		        +gamma,
		        -gamma,
		        -gamma,
		        +gamma,
		        +gamma,
		        -gamma,
		        +gamma,
		        -gamma,
		]
		bb = [
		        +gamma,
		        -gamma,
		        +gamma,
		        -gamma,
		        -gamma,
		        +gamma,
		        +gamma,
		        -gamma,
		]
		cc = [
		        +gamma,
		        -gamma,
		        +gamma,
		        -gamma,
		        +gamma,
		        -gamma,
		        -gamma,
		        +gamma,
		]

		# icosahedron: rest 12 faces
		aa.extend([+zeta, -zeta, +zeta, -zeta, 0, 0, 0, 0, +eta, -eta, +eta, -eta])
		bb.extend([+eta, -eta, -eta, +eta, +zeta, -zeta, +zeta, -zeta, 0, 0, 0, 0])
		cc.extend([0, 0, 0, 0, +eta, -eta, -eta, +eta, +zeta, -zeta, -zeta, +zeta])
	else:
		errors += 1
		logging.error('Invalid numFaces. Should be either: {4: tetrahedron, 6: cube, 8: octahedron, 12: dodecahedron, 20: icosahedron}')
		return (None)


#	if errors==0:
	if edge > 0.:
		ri = (edge / 2.) / tan(pi / p[numFaces]) * t
		rm = (edge / 2.) * cos(pi / p[numFaces]) / sin(pi / h[numFaces])
		rc = (edge / 2.) * tan(pi / q[numFaces]) * t
	elif ri > 0.:
		edge = 2 * ri * tan(pi / p[numFaces]) / t
		rm = (edge / 2.) * cos(pi / p[numFaces]) / sin(pi / h[numFaces])
		rc = (edge / 2.) * tan(pi / q[numFaces]) * t
	elif rm > 0.:
		edge = 2 * rm * sin(pi / h[numFaces]) / cos(pi / p[numFaces])
		ri = (edge / 2.) / tan(pi / p[numFaces]) * t
		rc = (edge / 2.) * tan(pi / q[numFaces]) * t
	elif rc > 0.:
		edge = 2 * rc / (tan(pi / q[numFaces]) * t)
		ri = (edge / 2.) / tan(pi / p[numFaces]) * t
		rm = (edge / 2.) * cos(pi / p[numFaces]) / sin(pi / h[numFaces])
	elif volume > 0.:
		edge = (24 * volume * (tan(pi / p[numFaces])**2) / (t * numFaces * p[numFaces]))**(1 / 3)
		ri = (edge / 2.) / tan(pi / p[numFaces]) * t
		rm = (edge / 2.) * cos(pi / p[numFaces]) / sin(pi / h[numFaces])
		rc = (edge / 2.) * tan(pi / q[numFaces]) * t

	dd = [ri] * numFaces
	if not r:
		r = min(dd) / 2.
	platonic = potentialblock(material=material, a=aa, b=bb, c=cc, d=array(dd) - r, r=r, R=R, mask=mask, isBoundary=isBoundary, fixed=fixed, color=color)

	platonic.shape.edge = edge
	platonic.shape.ri = ri  # Add inradius attribute from Python (visible only to Python objects/functions)
	platonic.shape.rm = rm  # Add midradius attribute from Python (visible only to Python objects/functions)
	platonic.shape.rc = rc  # Add circumradius attribute from Python (visible only to Python objects/functions)

	platonic.state.pos = center
	platonic.state.ori = platonic.shape.orientation
	return platonic
