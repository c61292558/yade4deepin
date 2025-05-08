# -*- encoding=utf-8 -*-
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
# 2023 © Robert Caulk <rob.caulk@gmail.com>

# Benchmarks from Chung and Ooi (2011)
# Test No 3: Normal contact with different restitution coefficients
# Units: SI (m, N, Pa, kg, sec)

# script execution:
# yade Test3.py (with default arguments, else:)
# yade Test3.py [aluminium_oxide_oxide|cast_iron] [en] [use_en]
# example: yade Test3.py aluminium_oxide 1.0 False (equivalent to default)

# ------------------------------------------------------------------------------------------
import math
from yade import plot

material = str(sys.argv[1]) if len(sys.argv) > 1 else 'aluminium_oxide'  # options are 'aluminium_oxide'|'cast_iron'
en = float(str(sys.argv[2])) if len(sys.argv) > 2 else 1.0  # options are 0.2, 0.4, 0.6, 0.8, 1.0
use_en = (sys.argv[3]) if len(sys.argv) > 3 else 'False'  # optional parameter, whether to define en/es or betan/betas
vn_ini = 3.9

try:
	os.mkdir('outputData')
except:
	pass  # will pass if folder already exists

# ------------------------------------------------------------------------------------------
# Material properties
if material == 'aluminium_oxide':
	spheresMat = O.materials.append(FrictMat(young=3.8e11, poisson=0.23, density=4000, frictionAngle=atan(0.0)))
	wallsMat = O.materials.append(
	        FrictMat(young=3.8e14, poisson=0.23, density=4000, frictionAngle=atan(0.0))
	)  # The wall/facet elements are 1000 stiffer than the sphere
elif material == 'cast_iron':
	spheresMat = O.materials.append(FrictMat(young=1.0e11, poisson=0.25, density=7000, frictionAngle=atan(0.0)))
	wallsMat = O.materials.append(
	        FrictMat(young=1.0e14, poisson=0.25, density=7000, frictionAngle=atan(0.0))
	)  # The wall/facet elements are 1000 stiffer than the sphere
else:
	print('The material should be either aluminium_oxide or cast_iron')

# ------------------------------------------------------------------------------------------
# Create boundary as wall/facet
r = 0.0025

O.bodies.append(sphere([-r, 0, 0], material=spheresMat, radius=r))
O.bodies.append(wall(position=0, axis=0, material=wallsMat))  # Use walls
#O.bodies.append(facet([[0,-3*r,-3*r],[0,-3*r,3*r],[0,3*r,0]],material=wallsMat,fixed=True,color=[1,0,0])) # Use facets

b0 = O.bodies[0]
b1 = O.bodies[1]

b0.state.vel[0] = vn_ini

hasHadContact = False  # Whether the two bodies have ever been in contact

if use_en == 'True':  # Here I use 'True' as a string
	ip2 = Ip2_FrictMat_FrictMat_MindlinPhys(en=en)
else:
	# Option 1: Traditional formula of beta-en
	#	betan=-log(en)/(pi**2+(log(en)**2))**0.5	# This formula is only correct if we allow viscous attraction.

	# Option 2: Thornton et al (2013) specifically for a Hertzian model with no attractive forces
	h1 = -6.918798
	h2 = -16.41105
	h3 = 146.8049
	h4 = -796.4559
	h5 = 2928.711
	h6 = -7206.864
	h7 = 11494.29
	h8 = -11342.18
	h9 = 6276.757
	h10 = -1489.915

	alpha = en * (h1 + en * (h2 + en * (h3 + en * (h4 + en * (h5 + en * (h6 + en * (h7 + en * (h8 + en * (h9 + en * h10)))))))))
	betan = 0 if en == 1.0 else sqrt(1 / (1 - ((1 + en)**2) * exp(alpha)) - 1)  # This is noted as gamma in Thornton et al (2013)

	betas = betan
	ip2 = Ip2_FrictMat_FrictMat_MindlinPhys(betan=betan, betas=betas)

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
                [ip2],
                [Law2_ScGeom_MindlinPhys_Mindlin(preventGranularRatcheting=False)]
        ),
        NewtonIntegrator(gravity=(0, 0, 0), damping=0),
        PyRunner(command='addPlotData()', iterPeriod=1)
]

O.dt = 0.01 * RayleighWaveTimeStep()


# ------------------------------------------------------------------------------------------
# Collect history of data which will be plotted
def addPlotData():
	global hasHadContact
	if O.interactions.countReal() > 0:
		ii = O.interactions.nth(0)
		ii.phys.kno = Kno
		ii.phys.kso = Kso
		hasHadContact = True
	else:
		if hasHadContact == True and b0.state.pos[0] < -r:
			vn_fin = b0.state.vel[0]

			flag = 'w' if en == 1.0 else 'a'  # Rewrite file for the first value en=0.2 or else append
			name = '_en_es' if use_en == 'True' else '_betan_betas'

			with open('outputData/Test3_' + material + name + '.txt', flag) as myfile:
				myfile.write(str(en))
				myfile.write(' ')
				myfile.write(str(-vn_fin / vn_ini))
				myfile.write("\n")
			O.pause()


O.saveTmp()

# ------------------------------------------------------------------------------------------
# Visualisation
try:
	from yade import qt
	v = qt.View()
	v.sceneRadius = 5 * r
except:
	pass

O.run()
