# -*- encoding=utf-8 -*-
# 2022 © Vasileios Angelidakis <vasileios.angelidakis@ncl.ac.uk>

# Duplicate script of cubePBscaled.py, simplified using the "potential_utils" module.

from yade import pack
from potential_utils import *
import math

recordVTK = False
# ----------------------------------------------------------------------------------------------------------------------------
if recordVTK:
	# Create vtk directory
	import os
	import errno
	try:
		os.mkdir('./vtk/')
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass

O.engines = [
        ForceResetter(),
        InsertionSortCollider([PotentialBlock2AABB()], verletDist=0.01),
        InteractionLoop(
                [Ig2_PB_PB_ScGeom(twoDimension=False, unitWidth2D=1.0, calContactArea=True)],
                [Ip2_FrictMat_FrictMat_KnKsPBPhys(kn_i=1e8, ks_i=1e7, Knormal=1e8, Kshear=1e7, useFaceProperties=False, viscousDamping=0.2)],
                [Law2_SCG_KnKsPBPhys_KnKsPBLaw(label='law', neverErase=False, allowViscousAttraction=True, traceEnergy=False)]
        ),
        NewtonIntegrator(damping=0.0, exactAsphericalRot=True, gravity=[0, -9.81, 0]),
]

powderDensity = 2000
meanSize = 1.0
distanceToCentre = meanSize / 2.
wallThickness = 0.5 * meanSize

O.materials.append(
        FrictMat(young=-1, poisson=-1, frictionAngle=radians(0.0), density=powderDensity, label='frictionless')
)  #The normal and shear stifness values are determined in the IPhys functor, thus the young, poisson parameters of the FrictMat are not used.
lengthOfBase = 9.0 * meanSize
heightOfBase = 14.0 * meanSize
sp = pack.SpherePack()
mn, mx = Vector3(-0.5 * (lengthOfBase - wallThickness), 0.5 * meanSize, -0.5 *
                 (lengthOfBase - wallThickness)), Vector3(0.5 * (lengthOfBase - wallThickness), 7.0 * heightOfBase, 0.5 * (lengthOfBase - wallThickness))
R = sqrt(3.0) * distanceToCentre
sp.makeCloud(mn, mx, R, 0, 100, False)

r = 0.01 * meanSize
for s in sp:
	color = Vector3(random.random(), random.random(), random.random())
	# PB_utils.cuboid
	b = cuboid(material=O.materials['frictionless'], edges=[1.0, 1.0, 1.0], r=r, R=0.0, center=s[0], mask=1, isBoundary=False, fixed=False, color=color)
	b.state.ori = Quaternion((random.random(), random.random(), random.random()), random.random())  #s[2]
	O.bodies.append(b)

#Bottom faces of the box
r = 0.1 * wallThickness
color = [0, 0.5, 1]

# Centroids of particles forming the bottom face of the box
center = [
        [0, 0, 0], [-lengthOfBase / 3.0, 0, 0], [lengthOfBase / 3.0, 0, 0], [0, 0, -lengthOfBase / 3.0], [-lengthOfBase / 3.0, 0, -lengthOfBase / 3.0],
        [lengthOfBase / 3.0, 0, -lengthOfBase / 3.0], [0, 0, lengthOfBase / 3.0], [-lengthOfBase / 3.0, 0, lengthOfBase / 3.0],
        [lengthOfBase / 3.0, 0, lengthOfBase / 3.0]
]

edges = [lengthOfBase / 3., wallThickness, lengthOfBase / 3.]

lidID = len(O.bodies)
for i in range(0, 9):
	O.bodies.append(
	        cuboid(material=O.materials['frictionless'], edges=edges, r=r, R=0.0, center=center[i], mask=3, isBoundary=True, fixed=True, color=color)
	)

# Vertical faces A-B-C-D of the box
center = [
        [0, 0.5 * heightOfBase, 0.5 * lengthOfBase], [0, 0.5 * heightOfBase, -0.5 * lengthOfBase], [0.5 * lengthOfBase, 0.5 * heightOfBase, 0],
        [-0.5 * lengthOfBase, 0.5 * heightOfBase, 0]
]

edges = [
        [lengthOfBase, heightOfBase, wallThickness], [lengthOfBase, heightOfBase, wallThickness], [wallThickness, heightOfBase, lengthOfBase],
        [wallThickness, heightOfBase, lengthOfBase]
]

transparentID = len(O.bodies)
for i in range(0, 4):
	O.bodies.append(
	        cuboid(material=O.materials['frictionless'], edges=edges[i], r=r, R=0.0, center=center[i], mask=3, isBoundary=True, fixed=True, color=color)
	)
O.bodies[transparentID].shape.wire = True

escapeNo = 0


def myAddPlotData():
	global escapeNo
	global wallThickness
	global meanSize
	uf = utils.unbalancedForce()
	if isnan(uf):
		uf = 1.0
	KE = utils.kineticEnergy()

	for b in O.bodies:
		if b.state.pos[1] < -5.0 * meanSize and len(b.state.blockedDOFs) == 0:  #i.e. fixed==False
			escapeNo = escapeNo + 1
			O.bodies.erase(b.id)
	if O.iter > 15000:
		removeLid()
	plot.addData(timeStep1=O.iter, timeStep2=O.iter, timeStep3=O.iter, timeStep4=O.iter, time=O.time, unbalancedForce=uf, kineticEn=KE, outsideNo=escapeNo)


from yade import plot

plot.plots = {'timeStep1': ('unbalancedForce'), 'timeStep2': ('kineticEn'), 'time': ('outsideNo')}
#plot.plot() #Uncomment to view plots
O.engines = O.engines + [PyRunner(iterPeriod=10, command='myAddPlotData()')]


def removeLid():
	global lidID
	if (O.bodies[lidID]):
		O.bodies.erase(lidID)


O.dt = 0.2 * sqrt(0.3 * O.bodies[0].state.mass / 1.0e8)

# ----------------------------------------------------------------------------------------------------------------------------
if recordVTK:
	# Take snapshots
	#	O.engines=O.engines+[SnapshotEngine(iterPeriod=100,fileBase='vtk/jenga-',firstIterRun=10,label='snapshooter')]

	# Export VTK results
	from yade import export
	vtkExporter = export.VTKExporter('./vtk/cubePBscaled')

	def vtkExport():
		vtkExporter.exportPotentialBlocks(what=dict(n='b.id'))

	O.engines = O.engines + [PyRunner(iterPeriod=500, command='vtkExport()')]

from yade import qt

qt.Controller()
v = qt.View()

O.saveTmp()
