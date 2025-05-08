# -*- encoding=utf-8 -*-
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
# 2023 © Robert Caulk <rob.caulk@gmail.com>

# Benchmarks from Chung and Ooi (2011)
# Test No 1: Elastic normal impact of two identical spheres
# Units: SI (m, N, Pa, kg, sec)

# script execution:
# yade Test1.py (with default arguments, else:)
# yade Test1.py [glass|limestone]
# example: yade Test1.py glass (equivalent to default)

# ------------------------------------------------------------------------------------------
from yade import plot

material = str(sys.argv[1]) if len(sys.argv) > 1 else 'glass'  #options are 'glass'|'limestone'

try:
	os.mkdir('outputData')
except:
	pass  # will pass if folder already exists

# ------------------------------------------------------------------------------------------
# Material properties
if material == 'glass':
	spheresMat = O.materials.append(FrictMat(young=4.8e10, poisson=0.20, density=2800, frictionAngle=atan(0.35)))
elif material == 'limestone':
	spheresMat = O.materials.append(FrictMat(young=2.0e10, poisson=0.25, density=2500, frictionAngle=atan(0.35)))
else:
	print('The material should be either glass or limestone')

# ------------------------------------------------------------------------------------------

r = 0.010
O.bodies.append(sphere([-r, 0, 0], material=spheresMat, radius=r))
O.bodies.append(sphere([r, 0, 0], material=spheresMat, radius=r))

b0 = O.bodies[0]
b1 = O.bodies[1]

b0.state.vel[0] = 10
b1.state.vel[0] = -10

hasHadContact = False  # Whether the two bodies have ever been in contact

# ------------------------------------------------------------------------------------------
# Simulation loop
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_MindlinPhys()], [Law2_ScGeom_MindlinPhys_Mindlin()]),
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
		fn = i.phys.normalForce.norm()
		un = i.geom.penetrationDepth
		hasHadContact = True
	else:
		if hasHadContact == True and b0.state.pos[0] < -r:
			plot.saveDataTxt('outputData/Test1_' + material + '.txt')
			plot.plot(noShow=True).savefig('outputData/Test1_' + material + '.png')
			O.pause()
	plot.addData(Dn=un * 10**6, Fn=fn / 1000)  #i=O.iter, t=O.time, # Write the output data in μm and kN


plot.plots = {'Dn': ('Fn')}  #, 't':('Fn')
O.saveTmp()

addPlotData()

try:  #if 'QT5' in features
	from yade import qt
	v = qt.View()
except:
	pass

O.run()
