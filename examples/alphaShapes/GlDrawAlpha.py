# Demonstrates alpha shapes, change the parameters {alpha, shrinkedAlpha, fixedAlpha} and visualize the result interactively

# First, create a non-trivial shape: build a cube, make a spherical cut through one corner, and dig a tunnel into it
# _________________________________________________________________________________________________________________

import random

N = 30
for i in range(N):
	for j in range(N):
		for k in range(N):
			O.bodies.append(
			        sphere((i + 0.01 * random.random(), j + 0.01 * random.random(), k + 0.01 * random.random()), 0.5, color=(0.7, 0.7, 0.1))
			)

n = Vector3(1, 1, 1)
n.normalize()
for b in O.bodies:
	if (b.state.pos).norm() < N / 2:
		O.bodies.erase(b.id)
	if (b.state.pos - b.state.pos.dot(n) * n).norm() < 2.5:
		O.bodies.erase(b.id)

# Now detect the contour with an alpha shape, the next two methods are equivalent
# _______________________________________________________________________________

from yade import qt

alpha = 4
shrinkedAlpha = (sqrt(alpha) - 0.5)**2

if True:  #method 1, the functor will instantiate its own TesselationWrapper
	graph = GlExtra_AlphaGraph()
	graph.alpha = alpha
	graph.shrinkedAlpha = shrinkedAlpha
	graph.fixedAlpha = False
	qt.Renderer().extraDrawers = [graph]
	# later, change the parameters with "graph.alpha= somethingElse"

else:  #method 2, plug an existing TW into the Gl functor
	TW = TesselationWrapper()
	#TW.triangulate()
	TW.getAlphaGraph(alpha=alpha, shrinkedAlpha=shrinkedAlpha, fixedAlpha=False)
	qt.Renderer().extraDrawers = [GlExtra_AlphaGraph(tesselationWrapper=TW)]
	# later, change the parameters with "TW.getAlphaGraph(alpha=other, ...)"

v = qt.View()
