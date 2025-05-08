# -*- encoding=utf-8 -*-
# 2022 © Vasileios Angelidakis <vasileios.angelidakis@ncl.ac.uk>

# Gravity deposition of regular polyhedra (Platonic solids) in a cylindrical packing, using the Potential Blocks

from yade import pack
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
# Create cylindrical packing of non-intersecting platonic solids with the same circumradius
rc = 0.05
radius = 0.5  #radius of the cylinder
height = 1.0  #height  of the cylinder

sp = pack.SpherePack()
mn = -Vector3(radius, radius, 0)
mx = Vector3(radius, radius, height)
sp.makeCloud(mn, mx, rc, 0, 200, False)

cyl = pack.inCylinder((0, 0, 0), (0, 0, 0.999 * height), 0.999 * radius)  # Create predicate inCylinder((bottomPoint,topPoint),radius)
sp = pack.filterSpherePack(cyl, sp, True)  # Filter initial cloud to retain only the bodies inside the predicate

r = 0.3 * rc
for s in sp:
	pb = platonic_solid(
	        material=m, numFaces=random.choice([4, 6, 8, 12]), rc=rc, r=r, R=0.0, center=s[0], mask=1, isBoundary=False, fixed=False
	)  # FIXME: Icosahedra are excluded for now
	pb.state.ori = Quaternion((random.random(), random.random(), random.random()), random.random())  #s[2]
	O.bodies.append(pb)

plates = cylindricalPlates(
        material=m, radius=0.5, height=1.0, thickness=0.05, numFaces=12, r=0.005, R=0.0, mask=3, isBoundary=False, fixed=True, lid=[False, True]
)
O.bodies.append(plates)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Visualisation
try:
	from yade import qt
	qt.Controller()
	v = qt.View()
	v.ortho = True
	v.sceneRadius = 20

	v.eyePosition = Vector3(0, -2, 0.5)
	v.viewDir = Vector3(0, 1, 0)
	v.upVector = Vector3(0, 0, 1)

	rndr = yade.qt.Renderer()
except:
	pass

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Timestep
O.dt = 2 * sqrt(O.bodies[0].state.mass / Kn)

O.saveTmp()
O.run()
