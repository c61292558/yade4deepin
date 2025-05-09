# -*- encoding=utf-8 -*-

cylHt, cylRd = 1, .2
nSpheres = 2e4


def unitSquare():
	"""Return square composed of 2 facets in the xy plane, centered at zero with unit extents."""
	import gts
	vv = [gts.Vertex(1, 1, 0), gts.Vertex(-1, 1, 0), gts.Vertex(-1, -1, 0), gts.Vertex(1, -1, 0)]
	ee = [gts.Edge(vv[0], vv[1]), gts.Edge(vv[1], vv[2]), gts.Edge(vv[2], vv[3]), gts.Edge(vv[3], vv[0]), gts.Edge(vv[0], vv[2])]
	surf = gts.Surface()
	surf.add(gts.Face(ee[0], ee[1], ee[4]))
	surf.add(gts.Face(ee[2], ee[3], ee[4]))
	return surf


def unitCylinder(nDiv=24):
	"""Returns GTS surface approximating cylinder (without caps), with
	of height 2 and radius 2, centered at origin, axis coincident with
	the z-axis.

	:param int nDiv: polyhedron approximating circle.
	"""
	import numpy
	from yade import pack
	thetas = numpy.linspace(0, 2 * pi, nDiv, endpoint=True)
	ptsBase = [Vector3(cos(th), sin(th), -1) for th in thetas]
	ptsTop = [p + Vector3(0, 0, 2) for p in ptsBase]
	return pack.sweptPolylines2gtsSurface([ptsBase, ptsTop])


from yade import pack, timing

cyl = unitCylinder()
sq = unitSquare()
sq.translate(0, 0, -1)
cyl.copy(sq)
cyl.scale(cylRd, cylRd, .5 * cylHt)
cyl.rotate(1, 0, 0, -pi / 4)  # 45° anti-colckwise in the yz plane
# calling gtsSurface2Facets with just "cyl" (without constructing the faces tuple) ignores 2 faces that were copy'd before; bug in pygts?
cylIds = O.bodies.append(pack.gtsSurface2Facets(cyl))
sp = pack.SpherePack()
wd = cylRd * sqrt(2)
rMean = (.2 * wd * wd * cylHt / (nSpheres * (4 / 3.) * pi))**(1 / 3.)
print('Generating cloud…')
sp.makeCloud((-wd / 2, -wd / 2, -.5 * cylHt), (wd / 2, wd / 2, .5 * cylHt), rMean, 0, int(nSpheres), False)
sp.rotate((1, 0, 0), -pi / 4)
O.bodies.append([sphere(s[0], s[1]) for s in sp])

O.engines = [
        ForceResetter(),
        InsertionSortCollider([
                Bo1_Sphere_Aabb(),
                Bo1_Facet_Aabb(),
        ]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        RotationEngine(rotateAroundZero=True, zeroPoint=(0, 0, 0), rotationAxis=(0, 1, 1), angularVelocity=30 * (2 * pi / 60), ids=cylIds, label='rotor'),
        NewtonIntegrator(damping=.3, gravity=(0, 0, -1e3)),  # gravity artificially high, to make it faster going ;-)
]
O.dt = PWaveTimeStep()
O.stopAtIter = int(2 * pi / (rotor.angularVelocity * O.dt))
O.timingEnabled = True
timing.reset()

from yade import qt

qt.Controller()
qt.View()
