# -*- encoding=utf-8 -*-
### Simpificated buldozer simulation with VTK recorder
### vtk-files are saved in /tmp directory with names buldozer-*.vtk
from numpy import linspace
from numpy import arange
import gts
import itertools
from yade import geom, pack

###Initial Data
numKnifeParts = 10
radiusKnife = 1
lengthKnife = 2
buldozerHeight = 1.2
radiusSph = 0.05
numBoxes = Vector3(8, 5, 2)
gapBetweenBoxes = 0.05
sizeBox = (lengthKnife - (numBoxes[1] - 1) * gapBetweenBoxes) / numBoxes[1]

## PhysicalParameters
Density = 2400
frictionAngle = radians(35)
tc = 0.001
en = 0.3
es = 0.3

## Materials
facetMat = O.materials.append(ViscElMat(frictionAngle=frictionAngle, tc=tc, en=en, et=es))
sphereMat = O.materials.append(ViscElMat(density=Density, frictionAngle=frictionAngle, tc=tc, en=en, et=es))

### Creating the Buldozer Knife
### from facets, using GTS
Knife = []
for i in linspace(pi, pi * 3 / 2, num=numKnifeParts, endpoint=True):
	Knife.append(Vector3(radiusKnife * cos(i), 0, radiusKnife * sin(i)))

KnifeP = [Knife, [p + Vector3(0, lengthKnife, 0) for p in Knife]]
KnifePoly = pack.sweptPolylines2gtsSurface(KnifeP, threshold=1e-4)
KnifeIDs = O.bodies.append(pack.gtsSurface2Facets(KnifePoly, color=(1, 0, 0), wire=False, material=facetMat))

KnifeIDs += O.bodies.append(
        geom.facetBox(
                (-lengthKnife / 2 - radiusKnife, lengthKnife / 2, -radiusKnife + buldozerHeight / 2), (lengthKnife / 2, lengthKnife / 2, buldozerHeight / 2.),
                wallMask=47,
                color=(0, 1, 0),
                wire=False
        )
)

KnifeIDs += O.bodies.append(
        geom.facetBox(
                (-lengthKnife / 2 - radiusKnife - lengthKnife / 4., lengthKnife / 2, -radiusKnife + buldozerHeight * 3. / 2. - buldozerHeight / 4.),
                (lengthKnife / 4., lengthKnife / 3., buldozerHeight / 4.),
                wallMask=47,
                color=(0, 0, 1),
                wire=False
        )
)

O.bodies.append(
        geom.facetBox((0, 0, radiusKnife), (lengthKnife * 3, lengthKnife * 3, lengthKnife), wallMask=16, color=(1, 1, 1), wire=False, material=facetMat)
)

### Creating the material for buldozer
colorsph1 = Vector3(120, 234, 150)
colorsph2 = Vector3(0, 0, 1)

colorsph1.normalize()
colorsph2.normalize()
colorSph = colorsph1
for xyz in itertools.product(arange(0, numBoxes[0]), arange(0, numBoxes[1]), arange(0, numBoxes[2])):
	ids_spheres = O.bodies.appendClumped(
	        pack.regularHexa(
	                pack.inEllipsoid(
	                        (
	                                xyz[0] * (sizeBox + gapBetweenBoxes), xyz[1] * (sizeBox + gapBetweenBoxes) + sizeBox * 0.5, xyz[2] *
	                                (sizeBox + gapBetweenBoxes) - radiusKnife + sizeBox * 0.6
	                        ), (sizeBox / 2, sizeBox / 2, sizeBox / 2)
	                ),
	                radius=radiusSph,
	                gap=0,
	                color=colorSph,
	                material=sphereMat
	        )
	)
	if (colorSph == colorsph1):
		colorSph = colorsph2
	else:
		colorSph = colorsph1

O.dt = .2 * tc

O.engines = [
        ForceResetter(),
        InsertionSortCollider([
                Bo1_Sphere_Aabb(),
                Bo1_Facet_Aabb(),
        ]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
                [Law2_ScGeom_ViscElPhys_Basic()],
        ),
        TranslationEngine(translationAxis=[1, 0, 0], velocity=2, ids=KnifeIDs),  # Buldozer motion
        NewtonIntegrator(damping=0, gravity=[0, 0, -9.8]),
        VTKRecorder(iterPeriod=1000, fileName='/tmp/bulldozer-', recorders=['spheres', 'facets'])
]

O.saveTmp()
from yade import qt

qt.Controller()
qt.View()
r = qt.Renderer()
r.lightPos = Vector3(0, 0, 50)
O.run(20000)
#qt.makeSimulationVideo('/tmp/buldozer.ogg',iterPeriod=1000,fps=30)
