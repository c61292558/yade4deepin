# -*- encoding=utf-8 -*-
# 2022 © Vasileios Angelidakis <vasileios.angelidakis@ncl.ac.uk>

# Model multi-drum column under harmonic excitation, using the Potential Blocks

from yade import pack
from potential_utils import *
import math, os, errno
from numpy import array

# -----------------------------------------------------------------------------
m = O.materials.append(FrictMat(young=-1, poisson=-1, frictionAngle=atan(0.4), density=2000, label='frictionless'))

Kn = 1e8  #Pa/m
Ks = Kn * 2 / 3  #Pa/m

# -----------------------------------------------------------------------------
# Engines
O.engines = [
        ForceResetter(),
        InsertionSortCollider([PotentialBlock2AABB()], verletDist=0.01, avoidSelfInteractionMask=2),
        InteractionLoop(
                [Ig2_PB_PB_ScGeom(twoDimension=False, unitWidth2D=1.0, calContactArea=True)],
                [Ip2_FrictMat_FrictMat_KnKsPBPhys(kn_i=Kn, ks_i=Ks, Knormal=Kn, Kshear=Ks, useFaceProperties=False, viscousDamping=0.3)],
                [Law2_SCG_KnKsPBPhys_KnKsPBLaw(label='law', neverErase=False, allowViscousAttraction=False)]
        ),
        NewtonIntegrator(damping=0.0, exactAsphericalRot=True, gravity=[0, 0, -9.81], label='newton'),
]

# -----------------------------------------------------------------------------
# Create drum column with variable cross-sections

height = 10.08
maxR = 1.65
minR = 1.28

#drumHeight=0.25
drumNo = 10
drumHeight = height / drumNo

for i in range(0, drumNo):
	zBtm = (1 - (i + 0) * drumHeight / height)
	zTop = (1 - (i + 1) * drumHeight / height)
	btmRadius = minR + zBtm * (maxR - minR)
	topRadius = minR + zTop * (maxR - minR)
	b = prism(
	        material=m,
	        radius1=btmRadius,
	        radius2=topRadius,
	        thickness=drumHeight,
	        numFaces=10,
	        r=0.01,
	        R=0.0,
	        center=None,
	        color=[-1, -1, -1],
	        mask=1,
	        isBoundary=False,
	        fixed=False
	)
	b.state.pos = b.shape.position + [0, 0, i * (drumHeight)] + [0, 0, drumHeight / 2]
	O.bodies.append(b)

maxR = 1.28
minR = 1.65
for i in range(drumNo, drumNo + 2):
	zBtm = (1 - (i + 1) * drumHeight / height)
	zTop = (1 - (i + 0) * drumHeight / height)
	btmRadius = maxR + zBtm * (maxR - minR)
	topRadius = maxR + zTop * (maxR - minR)
	b = prism(
	        material=m,
	        radius1=topRadius,
	        radius2=btmRadius,
	        thickness=drumHeight,
	        numFaces=10,
	        r=0.01,
	        R=0.0,
	        center=None,
	        color=[-1, -1, -1],
	        mask=1,
	        isBoundary=False,
	        fixed=False
	)
	b.state.pos = b.shape.position + [0, 0, i * (drumHeight)] + [0, 0, drumHeight / 2]
	O.bodies.append(b)

# Top plate
length = 3 * minR
thickness = maxR / 1.5
pos = O.bodies[-1].state.pos + [0, 0, drumHeight / 2 + thickness / 2.]
pb = cuboid(material=m, edges=[length, length, thickness], r=0.01, R=0.0, center=pos, mask=1, isBoundary=False, fixed=False, color=[1, 0, 0])
O.bodies.append(pb)

# Bottom plate
length = 10 * minR
thickness = maxR
pos = [0, 0, -thickness / 2]
pb = cuboid(material=m, edges=[length, length, thickness], r=0.01, R=0.0, center=pos, mask=1, isBoundary=False, fixed=True, color=[0, 1, 0])
pb_ID = O.bodies.append(pb)

O.engines = O.engines + [HarmonicMotionEngine(A=[0.1, 0, 0], f=[3, 0, 0], fi=[pi / 2, 0, 0], ids=[pb_ID], label='hME')]  # Harmonic Motion Engine

# -----------------------------------------------------------------------------
# Visualisation
from yade import qt

qt.Controller()
v = qt.View()
v.sceneRadius = 20
v.ortho = True

v.eyePosition = Vector3(0, -23, 6)
v.upVector = Vector3(0, 0, 1)
v.viewDir = Vector3(0, 1, 0)

rndr = yade.qt.Renderer()
rndr.wire = True

# -----------------------------------------------------------------------------
# Time step
O.dt = 0.001 * sqrt(O.bodies[0].state.mass / Kn)

O.saveTmp()
O.run()
