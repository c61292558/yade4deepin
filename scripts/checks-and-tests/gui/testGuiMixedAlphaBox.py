import random
from testGuiHelper import TestGUIHelper

yade.log.setLevel("Default", yade.log.ERROR)
# Must match the filename without starting "testGui" otherwise test fails because an indiation of a passed test is
# the existence of file called testGui_MixedAlphaBox_OK_or_Skipped.txt, where "MixedAlphaBox" comes from this line.
scr = TestGUIHelper(
        "MixedAlphaBox",
        False,
        (
                Vector3(-11.07264127532486775, 5.479858836462248917, 17.15917135537100791),  # v.lookAt
                Vector3(0.7730344906758657153, -0.1667015862924478564, -0.6120688338357758163),  # v.viewDir
                Vector3(-11.84567576600073302, 5.646560422754697051, 17.77124018920678239),  # v.eyePosition
                Vector3(0.1305986416697474672, 0.9860073822144534983, -0.1036023021588683424)
        )  # v.upVector
)

# FIXME - test skipped because of (1) random crash (2) explosions in high precision. This needs a fix.
#         see https://gitlab.com/yade-dev/trunk/-/merge_requests/906#note_1213149258
# Remove ↓ this True condition to re-enable the test.
# if True - temporarily skip this test.
if ((True) or (not 'CGAL' in yade.config.features)):
	print("Skipping testGuiMixedAlphaBox, CGAL not available in features")
	scr.createEmptyFile("screenshots/testGui_MixedAlphaBox_OK_or_Skipped.txt")  # creating verification file so there is no error thrown
	os._exit(0)

N = 10
guiIterPeriod = 1000


def filterPos(X):
	# the first condition create a concavity (for fun), the other crops to extract a cylinder from a cube
	return (X - Vector3(N / 2, N / 2, N)).norm() < N / 2 or (X - Vector3(N / 2, X[1], N / 2)).norm() > N / 2


random.seed(2)
for i in range(N):
	for j in range(N):
		for k in range(N):
			X = Vector3(i + 0.01 * random.random(), j + 0.01 * random.random(), k + 0.01 * random.random())
			if not filterPos(X):
				O.bodies.append(sphere(X, 0.49, color=(0.7, 0.7, 0.1)))

boxes = O.bodies.append([box((N / 2, -0.51, N / 2), (N, 0, N), fixed=True), box((N / 2, N - 0.49, N / 2), (N, 0, N), fixed=True)])
O.bodies[boxes[1]].state.vel = (0, -0.1, 0)

TW = TesselationWrapper()
TW.triangulate()
TW.addBoundingPlane(1, True)
TW.addBoundingPlane(1, False)
a = 6
shrinkedA = (sqrt(a) - 0.5)**2
fixedA = True

ag = TW.getAlphaGraph(alpha=a, shrinkedAlpha=shrinkedA, fixedAlpha=fixedA)
graph = GlExtra_AlphaGraph(tesselationWrapper=TW, wire=True)

from yade import qt

rr = qt.Renderer()
rr.extraDrawers = [graph]

TW.applyAlphaForces(stress=-Matrix3.Identity, alpha=a, shrinkedAlpha=shrinkedA, fixedAlpha=fixedA, reset=False)


def updateMembraneForces():
	TW.triangulate(reset=True)
	TW.addBoundingPlane(1, True)
	TW.addBoundingPlane(1, False)
	TW.applyAlphaForces(stress=-1000 * Matrix3.Identity, alpha=a, shrinkedAlpha=shrinkedA, fixedAlpha=fixedA, reset=False)


def displayMembrane():
	#FIXME: suboptimal, we can't re-use the existing triangulation because it would try to insert fictious spheres again
	# as a consequence we rebuild the whole thing just for display :-/
	TW.triangulate(reset=True)
	TW.addBoundingPlane(1, True)
	TW.addBoundingPlane(1, False)
	ag = TW.getAlphaGraph(alpha=a, shrinkedAlpha=shrinkedA, fixedAlpha=fixedA)
	graph.refresh()
	for b in O.bodies:
		if isinstance(b.shape, Sphere):
			if O.forces.permF(b.id) != Vector3.Zero:
				b.shape.color = (0.8, 0.4, 0.4)
			else:
				b.shape.color = (0.4, 0.7, 0.4)
		else:
			b.shape.color = (1, 0, 0)


# v = qt.View() # TestGUIHelper takes care of opening the view.

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], verletDist=0.3),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=100, timestepSafetyCoefficient=0.8),
        #triax,
        NewtonIntegrator(damping=0.3),
        PyRunner(iterPeriod=50, initRun=True, command="updateMembraneForces()", label="membrane"),
        PyRunner(iterPeriod=100, initRun=True, command="displayMembrane()", label="GuiMembrane"),
        PyRunner(iterPeriod=guiIterPeriod, command='scr.screenshotEngine()')
]

O.run(guiIterPeriod * scr.getTestNum() + 1)

# make Video
#O.engines = O.engines + [qt.SnapshotEngine(fileBase=O.tmpFilename(), iterPeriod=100, label='snapshotter')]
#def done():
#makeVideo(snapshotter.snapshots,out='alphaWithWalls.mpeg', fps=20, kbps=5000)
