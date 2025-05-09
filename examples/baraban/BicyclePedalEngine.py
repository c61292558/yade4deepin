# -*- encoding=utf-8 -*-
"""
Engine applying the linear motion of bicycle pedal e.g. moving points around the axis without rotation
"""

import time

## PhysicalParameters
Density = 2400
frictionAngle = radians(35)
tc = 0.001
en = 0.3
es = 0.3

## Import wall's geometry
facetMat = O.materials.append(ViscElMat(frictionAngle=frictionAngle, tc=tc, en=en, et=es))  # **params sets kn, cn, ks, cs
sphereMat = O.materials.append(ViscElMat(density=Density, frictionAngle=frictionAngle, tc=tc, en=en, et=es))
from yade import ymport

fctIds = O.bodies.append(ymport.stl('baraban.stl', color=(1, 0, 0), material=facetMat))
## Spheres
sphereRadius = 0.2
nbSpheres = (10, 10, 10)
#nbSpheres = (50,50,50)
for i in range(nbSpheres[0]):
	for j in range(nbSpheres[1]):
		for k in range(nbSpheres[2]):
			x = (i * 2 - nbSpheres[0]) * sphereRadius * 1.1
			y = (j * 2 - nbSpheres[1]) * sphereRadius * 1.1
			z = (k * 2 - nbSpheres[2]) * sphereRadius * 1.1
			s = sphere([x, y, z], sphereRadius, material=sphereMat)
			O.bodies.append(s)

## Timestep
O.dt = .2 * tc

## Engines
O.engines = [
        ## Resets forces and momenta the act on bodies
        ForceResetter(),
        ## Using bounding boxes find possible body collisions.
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        ## Interactions
        InteractionLoop(
                ## Create geometry information about each potential collision.
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                ## Create physical information about the interaction.
                [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
                ## Constitutive law
                [Law2_ScGeom_ViscElPhys_Basic()],
        ),
        ## Apply gravity
        ## Cundall damping must been disabled!
        NewtonIntegrator(damping=0, gravity=[0, -9.81, 0]),
        ## Saving results
        #VTKRecorder(virtPeriod=0.04,fileName='/tmp/stlimp-',recorders=['spheres','facets']),
        ## Apply kinematics to walls
        BicyclePedalEngine(ids=fctIds, rotationAxis=[0, 0, 1], radius=1.2, angularVelocity=4.0)
]

from yade import qt

qt.View()
#O.saveTmp()
#O.run()
