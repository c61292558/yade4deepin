if ('CGAL' in features):
	import random
	random.seed(1)
	N = 30
	for i in range(N):
		for j in range(N):
			for k in range(N):
				O.bodies.append(
				        sphere(
				                (i + 0.01 * random.random(), j + 0.01 * random.random(), k + 0.01 * random.random()),
				                0.5,
				                color=(0.7, 0.7, 0.1)
				        )
				)
	if False:
		O.bodies.append(sphere((-1000.5, N / 2, N / 2), 1000, color=(0.7, 0.7, 0.1)))

	n = Vector3(1, 1, 1)
	n.normalize()
	for b in O.bodies:
		if (b.state.pos).norm() < N / 2:
			O.bodies.erase(b.id)
		elif (b.state.pos - b.state.pos.dot(n) * n).norm() < 2.5:
			O.bodies.erase(b.id)

	TW = TesselationWrapper()
	TW.triangulate()
	a = 4
	shrinkedA = (sqrt(a) - 0.5)**2

	ag = TW.getAlphaGraph(alpha=a, shrinkedAlpha=shrinkedA, fixedAlpha=False)

	if len(ag) != 63640:
		raise YadeCheckError("the number of segments changed")
# this one is not reproducible, random numbering of segments
#if ag[0]!=Vector3(16.50265245538970049,23.978674468408812,-0.496985303100643705):
#raise YadeCheckError("a vertex position changed", ag[0],"!=",Vector3(16.50265245538970049,23.978674468408812,-0.496985303100643705))

# for display, comment out
#graph = GlExtra_AlphaGraph(tesselationWrapper=TW,wire=False)
#graph.lighting=True
#from yade import qt
#rr = qt.Renderer().extraDrawers = [graph]
#v = qt.View()
else:
	print("CGAL not in features, skip")
