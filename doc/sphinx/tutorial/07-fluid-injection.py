# This script simulates the injection of a fluid in a localized region below immersed particles
# The simulation is periodic along z-axis.
# at first execution, the simulation starts by depositing the particles in the container then saves the scene before proceeding to injection
# further execution will reload the deposited layer and start injection directly to gain time
# WARNING: changes in some input parameters like dimensions of the box or number of particles will not be reflected as long as the saved state is present on disk,
#          remember to erase it to force a new generation, or set newSample=True below

from yade import pack, export
import yade.timing
from math import *
from pylab import rand
import os.path
import numpy as np
import matplotlib.pyplot as plt

O.periodic = True

# Dimensions of the box and of the injection zone
width = 0.8  #
height = 1
depth = 0.4
aperture = 0.05 * width

# number of spheres
numSpheres = 2000
# contact friction during deposition
compFricDegree = 10
# porosity of the initial cloud
porosity = 0.8
# Cundall's damping (zero recommanded)
damp = 0.0
# fluid viscosity
mu = 0.01
# flow rate at the inlet
injectedFlux = -0.001
# name of output folder
key = 'output0'

newSample = False  #turn this true if you want to generate new sample by pluviation each time you run the script

# Deduced mean size for generating the cloud and consistency check
filename = "init" + key + str(numSpheres)
volume = width * height * depth
meanRad = pow(volume * (1 - porosity) / (pi * (4 / 3.) * numSpheres), 1 / 3.)
if (meanRad * 6 > depth):
	print("INCOMPATIBLE SIZES. INCREASE DEPTH OR INCREASE NUM_SPHERES")

# if no deposited layer has been found (first execution), generate bodies
if not os.path.isfile(filename + ".yade") or newSample:  #we create new sample if it does not exist...
	O.cell.hSize = Matrix3(width, 0, 0, 0, 3. * height, 0, 0, 0, depth)
	O.materials.append(FrictMat(young=400000.0, poisson=0.5, frictionAngle=compFricDegree / 180.0 * pi, density=1600, label='spheres'))
	O.materials.append(FrictMat(young=400000.0, poisson=0.5, frictionAngle=radians(15), density=1000, label='walls'))
	lowBox = box(center=(width / 2.0, height, width / 2.0), extents=(width * 1000.0, 0, width * 1000.0), fixed=True, wire=False)
	O.bodies.append(lowBox)
	topBox = box(center=(width / 2.0, 2 * height + 4 * meanRad, width / 2.0), extents=(width * 1000.0, 0, width * 1000.0), fixed=True, wire=False)
	O.bodies.append(topBox)

	sp = pack.SpherePack()
	sp.makeCloud((0, height + 2 * meanRad, 0), (width, 2 * height + 2 * meanRad, depth), -1, .002, numSpheres, periodic=True, porosity=porosity, seed=2)
	O.bodies.append([sphere(s[0], s[1], color=(0.6 + 0.15 * rand(), 0.5 + 0.15 * rand(), 0.15 + 0.15 * rand()), material='spheres') for s in sp])
else:  #... else we re-use the previous one
	O.load(filename + ".yade")

# look better
Gl1_Sphere.stripes = True

# define the fluid solver with appropriate parameter (see PeriodicFlowEngine's documentation)
flow = PeriodicFlowEngine(
        dead=1,
        meshUpdateInterval=40,
        defTolerance=-1,
        permeabilityFactor=1.0,
        useSolver=3,
        duplicateThreshold=depth,
        wallIds=[-1, -1, 0, 1, -1, -1],
        bndCondIsPressure=[0, 0, 0, 1, 0, 0],
        viscosity=mu,
        label="flow"
)

newton = NewtonIntegrator(damping=damp, gravity=(0, -9.81, 0))

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], allowBiggerThanPeriod=1),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=1000, timestepSafetyCoefficient=0.3, defaultDt=utils.PWaveTimeStep(), label='timestepper'),
        flow,  #  <====== the solver is inserted here, for the moment it is "dead" (doing nothing)
        newton,
        # some recorders for post-processing
        PyRunner(command="flow.saveVtk(key)", iterPeriod=25, dead=1, label="vtkE"),
        VTKRecorder(recorders=["spheres", "velocity", "stress"], iterPeriod=25, dead=1, fileName=key + '/', label="vtkR")
]

########## if this is fresh execution, get static equilibrium and save the result for later use  #############
if not os.path.isfile(filename + ".yade") or newSample:
	O.run(1000, 1)
	while unbalancedForce() > 0.01:
		O.run(100, 1)
	# turn the recorders and the solver on
	vtkE.dead = vtkR.dead = flow.dead = 0
	#add a recorder and define what to plot
	O.engines = O.engines + [PyRunner(command=('myAddPlotData()'), iterPeriod=50, label="recorder")]
	O.save(filename + ".yade")

O.saveTmp()

########## define what to plot ##############
from yade import plot
import pylab

# First, find particle at the top of the sample (for evaluating initial height of the layer)
maxY = 0
for s in O.bodies:
	if isinstance(s.shape, Sphere):
		pos = s.state.pos
		if pos[1] > maxY:
			maxY = pos[1]


def myAddPlotData():
	index = flow.getCell(0.5 * width, height, 0.5 * depth)
	if index > 0:
		########## find particle at the top of the sample ##############
		simpleH = 0
		for s in O.bodies:
			if isinstance(s.shape, Sphere):
				pos = s.state.pos
				if pos[1] > simpleH:
					simpleH = pos[1]
		########## function to compute hf ##############
		cavityh = height
		for s in O.bodies:
			v = s.state.vel
			magvel = pow((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]), 0.5)
			if magvel > 0.14:
				pos = s.state.pos
				if pos[1] > cavityh:
					cavityh = pos[1]

		#########################
		plot.addData(
		        t=O.time,
		        i=O.iter,
		        p=flow.getPorePressure((0.5 * width, height, 0.5 * depth)),
		        q=-injectedFlux,
		        Ho=maxY - 1,
		        hf=cavityh - 1,
		        H=simpleH - 1
		)

		H = simpleH - 1
		hf = cavityh - 1


################ impose the costant flux  ###############


# In this function we find the elements of the mesh which have a face in the injection region, to distribute the inlet flux
# it is inserted in the solver below
def imposeFlux(valF):
	found = 0
	listF1 = []
	for k in range(flow.nCells()):
		if 0 in flow.getVertices(k) and flow.getCellCenter(k)[0] > ((width - aperture) / 2.) and flow.getCellCenter(k)[0] < ((width + aperture) / 2.):
			listF1.append(k)
	flow.clearImposedFlux()
	if len(listF1) == 0:
		flow.imposeFlux((0.5 * width, height, 0.5 * depth), valF)
	else:
		for k in listF1:
			flow.imposeFlux(flow.getCellBarycenter(k), valF / float(len(listF1)))


injectedFlux = -1  #this one is large and it will fluidize violently, see below for smoother evolutions

# very important: the flux needs to be imposed with this "hook", which plugs our custom function in the right place in the solving sequence
flow.blockHook = "imposeFlux(injectedFlux)"

###############  EXERCISE ###################
# 1- use trial and error in the shell (see the note below) to find approximately an upper-bound of "injectedFlux", below which there is only limited movements of the particles
# 2- implement in this script a progressive increase of the flux, starting from a small fraction of max value above, and exceeding it by a factor 2 at the end
# 3- choose an approriate plot to show that the pressure response is initially linear, then sublinear. Is the upper-bound from question 1 also an upper-bound for the linear response?
# 4- compare with fig. 10 from https://doi.org/10.1103/PhysRevE.94.052905 (also available at https://arxiv.org/pdf/1703.02319)
# 5- use paraview to generate a video similar to https://www.youtube.com/watch?v=gH585XaQEcY (it can be with a constant flux)

#NOTE:
# to change the injected flux interactively, change it in the global scope:
# globals()["injectedFlux"]=-0.03
# else we have two variables with the same name in different scopes and
