# -*- encoding=utf-8 -*-
# 2022 © Vasileios Angelidakis <vasileios.angelidakis@ncl.ac.uk>

# Showcase of shapes that can be generated using the potential_utils module

from yade import pack, utils
from potential_utils import *
import math, os, errno
from numpy import array

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Material Parameters
m = O.materials.append(FrictMat(young=-1, poisson=-1, frictionAngle=atan(0.5), density=2000, label='frictMat'))

Kn = 1e8  #Pa/m
Ks = Kn * 2 / 3  #Pa/m

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Engines
O.engines = [
        ForceResetter(),
        InsertionSortCollider([PotentialBlock2AABB()], verletDist=0.01, avoidSelfInteractionMask=2),
        InteractionLoop(
                [Ig2_PB_PB_ScGeom(twoDimension=False, unitWidth2D=1.0, calContactArea=True)],
                [Ip2_FrictMat_FrictMat_KnKsPBPhys(kn_i=Kn, ks_i=Ks, Knormal=Kn, Kshear=Ks, useFaceProperties=False, viscousDamping=0.2)],
                [Law2_SCG_KnKsPBPhys_KnKsPBLaw(label='law', neverErase=False, allowViscousAttraction=False)]
        ),
        NewtonIntegrator(damping=0.0, exactAsphericalRot=True, gravity=[0, 0, -9.81], label='newton'),
]

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Create Potential Blocks of various shapes
rc = 0.05
length = 20.0 * rc
r = 0.1 * rc
edges = [rc / sqrt(3), rc / sqrt(3), rc / sqrt(3)]

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.potentialblock: Generic Potential Block defined by a,b,c,d
pb = potentialblock(
        material=m,
        a=[1, -1, 0, 0, 0, 0],
        b=[0, 0, 1, -1, 0, 0],
        c=[0, 0, 0, 0, 1, -1],
        d=[edges[0] / 2. - r, edges[0] / 2. - r, edges[1] / 2. - r, edges[1] / 2. - r, edges[2] / 2. - r, edges[2] / 2. - r],
        r=r,
        R=0.0,
        mask=1,
        isBoundary=False,
        fixed=False,
        color=[1, 0, 0]
)  # Red cuboid
pb.state.pos = [-8 * rc, rc, 0]
O.bodies.append(pb)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.cuboid: Cuboidal Potential Block
pb = cuboid(material=m, edges=edges, r=r, R=0.0, center=[-8 * rc, 0, 0], mask=1, isBoundary=False, fixed=True, color=[0, 1, 0])  # Green cuboid
O.bodies.append(pb)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.platonic_solid: The five Platonic Solids
xPos = -7 * rc
for numFaces in [4, 6, 8, 12, 20]:
	pb = platonic_solid(material=m, numFaces=numFaces, rm=0.01, r=r, R=0.0, center=[xPos, 0.05 * length, 0], mask=1, isBoundary=False, fixed=False)
	O.bodies.append(pb)
	xPos = xPos + rc

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.platonic_solid: Dodecahedra defining their volume, edge, inradius, midradius and circumradius
numFaces = 12  #4, 6, 8, 12, 20
pb = platonic_solid(material=m, numFaces=numFaces, volume=0.00001, r=r, R=0.0, center=[-7 * rc, 0, 0], mask=1, isBoundary=False, fixed=False)
O.bodies.append(pb)  # Define volume
pb = platonic_solid(material=m, numFaces=numFaces, edge=0.01, r=r, R=0.0, center=[-6 * rc, 0, 0], mask=1, isBoundary=False, fixed=False)
O.bodies.append(pb)  # Define edge
pb = platonic_solid(material=m, numFaces=numFaces, ri=0.01, r=r, R=0.0, center=[-5 * rc, 0, 0], mask=1, isBoundary=False, fixed=False)
O.bodies.append(pb)  # Define inradius
pb = platonic_solid(material=m, numFaces=numFaces, rm=0.01, r=r, R=0.0, center=[-4 * rc, 0, 0], mask=1, isBoundary=False, fixed=False)
O.bodies.append(pb)  # Define inradius
pb = platonic_solid(material=m, numFaces=numFaces, rc=0.01, r=r, R=0.0, center=[-3 * rc, 0, 0], mask=1, isBoundary=False, fixed=False)
O.bodies.append(pb)  # Define circumradius

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.prism: Prisms with different numbers of sides
xPos = -8 * rc
for numFaces in [3, 4, 5, 6, 7, 8]:
	pb = prism(
	        material=m,
	        radius1=0.2 * rc,
	        radius2=0.2 * rc,
	        thickness=0.2 * rc,
	        numFaces=numFaces,
	        r=0.0,
	        R=0.0,
	        center=None,
	        color=[-1, -1, -1],
	        mask=1,
	        isBoundary=False,
	        fixed=False
	)
	pb.state.pos = [xPos, 2 * rc, 0]
	xPos = xPos + rc
	O.bodies.append(pb)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.prism: Cylinders with varying numbers of faces
xPos = -8 * rc
for numFaces in [15, 20, 25, 30, 35, 40]:
	pb = prism(
	        material=m,
	        radius1=0.2 * rc,
	        radius2=0.2 * rc,
	        thickness=0.3 * rc,
	        numFaces=numFaces,
	        r=0.0,
	        R=0.0,
	        center=None,
	        color=[-1, -1, -1],
	        mask=1,
	        isBoundary=False,
	        fixed=False
	)
	pb.state.pos = [xPos, 3 * rc, 0]
	xPos = xPos + rc
	O.bodies.append(pb)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# potential_utils.prism: Cylinders with varying radius
xPos = -8 * rc
for radius2 in [0.05 * rc, 0.1 * rc, 0.15 * rc, 0.2 * rc, 0.25 * rc, 0.3 * rc]:
	pb = prism(
	        material=m,
	        radius1=0.2 * rc,
	        radius2=radius2,
	        thickness=0.5 * rc,
	        numFaces=20,
	        r=0.0,
	        R=0.0,
	        center=None,
	        color=[-1, -1, -1],
	        mask=1,
	        isBoundary=False,
	        fixed=False
	)
	pb.state.pos = [xPos, 4 * rc, 0]
	xPos = xPos + rc
	O.bodies.append(pb)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Visualisation
try:
	from yade import qt
	qt.Controller()
	v = qt.View()
	v.ortho = False
	v.sceneRadius = 20

	v.eyePosition = Vector3(-0.25, 0.1, 0.67)
	v.viewDir = Vector3(0, 0, -1)
	v.upVector = Vector3(0, 1, 0)

	rndr = yade.qt.Renderer()
except:
	pass

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Timestep
O.dt = 2 * sqrt(O.bodies[0].state.mass / Kn)

O.saveTmp()
#O.run()
