# -*- encoding=utf-8 -*-
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
# 2023 © Robert Caulk <rob.caulk@gmail.com>

# Benchmarks from Chung and Ooi (2011)
# Test No 2: Elastic normal impact of a sphere with a rigid plane
# Units: SI (m, N, Pa, kg, sec)

# script execution:
# yade Test2.py (with default arguments, else:)
# yade Test2.py [aluminium_alloy_alloy|magnesium_alloy]
# example: yade Test2.py aluminium_alloy (equivalent to default)

# ------------------------------------------------------------------------------------------
from yade import plot

material = str(sys.argv[1]) if len(sys.argv) > 1 else 'aluminium_alloy'  #options are 'aluminium_alloy'|'magnesium_alloy'

try:
	os.mkdir('outputData')
except:
	pass  # will pass if folder already exists

# ------------------------------------------------------------------------------------------
# Material properties
if material == 'aluminium_alloy':
	spheresMat = O.materials.append(FrictMat(young=7.0e10, poisson=0.30, density=2699, frictionAngle=atan(0.0)))
	wallsMat = O.materials.append(
	        FrictMat(young=7.0e13, poisson=0.30, density=2699, frictionAngle=atan(0.0))
	)  # The wall/facet elements are 1000 stiffer than the sphere
elif material == 'magnesium_alloy':
	spheresMat = O.materials.append(FrictMat(young=4.0e10, poisson=0.35, density=1800, frictionAngle=atan(0.0)))
	wallsMat = O.materials.append(
	        FrictMat(young=4.0e13, poisson=0.35, density=1800, frictionAngle=atan(0.0))
	)  # The wall/facet elements are 1000 stiffer than the sphere
else:
	print('The material should be either aluminium_alloy or magnesium_alloy')

# ------------------------------------------------------------------------------------------
# Create boundary as wall/facet
r = 0.10

O.bodies.append(sphere([-r, 0, 0], material=spheresMat, radius=r))
O.bodies.append(wall(position=0, axis=0, material=wallsMat))  # Use walls
#O.bodies.append(facet([[0,-3*r,-3*r],[0,-3*r,3*r],[0,3*r,0]],material=wallsMat,fixed=True,color=[1,0,0])) # Use facets

b0 = O.bodies[0]
b1 = O.bodies[1]

b0.state.vel[0] = 0.2

hasHadContact = False  # Whether the two bodies have ever been in contact

# ------------------------------------------------------------------------------------------
# Calculate elastic contact parameters
Ea = O.materials[0].young
Va = O.materials[0].poisson
Ga = Ea / (2 * (1 + Va))

Eb = O.materials[1].young
Vb = O.materials[1].poisson
Gb = Eb / (2 * (1 + Vb))

E = Ea * Eb / ((1. - Va**2) * Eb + (1. - Vb**2) * Ea)
G = 1.0 / ((2 - Va) / Ga + (2 - Vb) / Gb)

Kno = 4.0 / 3.0 * E * sqrt(r)
Kso = 8.0 * sqrt(r) * G

# ------------------------------------------------------------------------------------------
# Simulation loop
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(), Bo1_Wall_Aabb()]),
        InteractionLoop(
                # handle facet+sphere and wall+sphere collisions
                [Ig2_Facet_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_MindlinPhys()],
                [Law2_ScGeom_MindlinPhys_Mindlin(preventGranularRatcheting=False)]
        ),
        NewtonIntegrator(gravity=(0, 0, 0), damping=0),
        PyRunner(command='addPlotData()', iterPeriod=1)
]

O.dt = 0.01 * RayleighWaveTimeStep()


# collect history of data which will be plotted
def addPlotData():
	global hasHadContact
	fn = 0
	un = 0
	if O.interactions.countReal() > 0:
		i = O.interactions.all()[0]
		i.phys.kno = Kno
		i.phys.kso = Kso
		fn = i.phys.normalForce.norm()
		un = i.geom.penetrationDepth
		hasHadContact = True
	else:
		if hasHadContact == True and b0.state.pos[0] < -r:
			plot.saveDataTxt('outputData/Test2_' + material + '.txt')
			plot.plot(noShow=True).savefig('outputData/Test2_' + material + '.png')
			O.pause()
	plot.addData(Dn=un * 10**6, Fn=fn / 1000)  #i=O.iter, t=O.time,


plot.plots = {'Dn': ('Fn')}  #, 't':('Fn')
O.saveTmp()

addPlotData()

try:  #if 'QT5' in features:
	from yade import qt
	v = qt.View()
except:
	pass

O.run()
