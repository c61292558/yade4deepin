# -*- encoding=utf-8 -*-

from yade import ymport

## Omega
o = Omega()

## PhysicalParameters
Density = 2400
frictionAngle = radians(35)
sphereRadius = 0.05
tc = 0.001
en = 0.3
es = 0.3

## Import wall's geometry
facetMat = O.materials.append(ViscElMat(frictionAngle=frictionAngle, tc=tc, en=en, et=es))
sphereMat = O.materials.append(ViscElMat(density=Density, frictionAngle=frictionAngle, tc=tc, en=en, et=es))

walls = O.bodies.append(ymport.stl('ring.stl', material=facetMat))


def fill_cylinder_with_spheres(sphereRadius, cylinderRadius, cylinderHeight, cylinderOrigin, cylinderSlope):
	spheresCount = 0
	for h in range(0, int(cylinderHeight / sphereRadius / 2)):
		for r in range(1, int(cylinderRadius / sphereRadius / 2)):
			dfi = asin(0.5 / r) * 2
			for a in range(0, int(6.28 / dfi)):
				x = cylinderOrigin[0] + 2 * r * sphereRadius * cos(dfi * a)
				y = cylinderOrigin[1] + 2 * r * sphereRadius * sin(dfi * a)
				z = cylinderOrigin[2] + h * 2 * sphereRadius
				s = sphere(
				        [x, y * cos(cylinderSlope) + z * sin(cylinderSlope), z * cos(cylinderSlope) - y * sin(cylinderSlope)],
				        sphereRadius,
				        material=sphereMat
				)
				o.bodies.append(s)
				spheresCount += 1
	return spheresCount


# Spheres
spheresCount = 0
spheresCount += fill_cylinder_with_spheres(sphereRadius, 0.5, 0.10, [0, 0, 0], radians(0))
print("Number of spheres: %d" % spheresCount)

## Engines
o.engines = [
        ## Resets forces and momenta the act on bodies
        ForceResetter(),

        ## Using bounding boxes find possible body collisions.
        InsertionSortCollider([
                Bo1_Sphere_Aabb(),
                Bo1_Facet_Aabb(),
        ]),
        # Interactions
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
        ## Apply kinematics to walls
        ## angularVelocity = 0.73 rad/sec = 7 rpm
        RotationEngine(ids=walls, rotationAxis=[0, 0, 1], rotateAroundZero=True, angularVelocity=0.73)
]

for b in O.bodies:
	if isinstance(b.shape, Sphere):
		b.state.blockedDOFs = 'zXY'

O.dt = 0.02 * tc

O.saveTmp('init')

from yade import qt

renderer = qt.Renderer()
renderer.wire = True
#qt.Controller()
qt.View()
O.run()
