if ('CGAL' in features):
	import numpy as np

	# Test TesselationWrapper for finding volumes and deformations in sphere systems
	# Some of the steps below have no effect (inserting walls, inserting+deleting spheres) but they were producing errors in old versions
	# we make sure TW remains immune to them

	O.engines = [ForceResetter(), NewtonIntegrator()]

	# 1. Insert stuff that should be disregarded by TesselationWrapper
	O.bodies.append(wall(position=0, axis=0))

	# 2. Generate spheres
	import random
	random.seed(1)
	N = 6
	shear = 0.1
	shift = 0.01
	for i in range(N):
		for j in range(N):
			for k in range(N):
				O.bodies.append(
				        sphere(
				                (i + shift * random.random(), j + shift * random.random(), k + shift * random.random()),
				                0.5 * (1 - shift * (1 + random.random())),
				                color=(0.7, 0.7, 0.1)
				        )
				)
	# 3. Another non-sphere
	O.bodies.append(wall(position=N * 1.1, axis=1))

	# 4. Erase some bodies, invalidating some id's
	n = Vector3(1, 1, 1)
	n.normalize()
	nErased = 0
	for b in O.bodies:
		if isinstance(b.shape, Sphere) and (b.state.pos).norm() < N / 2:
			O.bodies.erase(b.id)
			nErased += 1

	# 5. Set motion and analyze
	for b in O.bodies:
		b.state.blockedDOFs = "xyzXYZ"
		b.state.vel = Vector3(shear, 0, 0) * b.state.pos[1]  # with dt=1, it should give a d(u_x)/dy = shear

	O.dynDt = False
	O.dt = 1

	TW = TesselationWrapper()
	TW.setState(0)
	O.step()
	TW.setState(1)
	vp = TW.getVolPoroDef(False)
	vpd = TW.getVolPoroDef(True)

	TW.defToVtk("checkNoCrash.vtk")  # we don't check content though...
	if not os.path.isfile("checkNoCrash.vtk"):
		raise YadeCheckError("no output file")
	os.remove("checkNoCrash.vtk")

	vv1 = vp["vol"]
	vv2 = vpd["vol"]
	vv3 = [TW.volume(b.id) for b in O.bodies]
	vv4 = [TW.volume(id) for id in range(len(O.bodies))]  # not the same as vv2 after erasing bodies
	deformations = vpd["def"]

	sums = [np.sum(xx) for xx in [vv1, vv2, vv3, vv4]]

	if (np.abs(sums[0] - sums[1]) > 1e-10 or np.abs(sums[0] - sums[2]) > 1e-10 or np.abs(sums[0] - sums[2]) > 1e-10):
		raise YadeCheckError("total volumes are different: " + str(sums))

	vTol = 1e-4 * N**3 + (6 * shift * (N + 2 * shift)**2)  # we have random fluctuations if shifts are large, add it to tolerance
	if (np.abs(sums[0] - (N**3)) > vTol):
		raise YadeCheckError("large error on volumes: " + str(sums[0]))

	F = Matrix3.Zero
	for f in deformations:
		F += f
	F /= N**3 - nErased
	if (np.abs(F[0, 1] - O.dt * shear)) > 1e-10:
		raise YadeCheckError("large error on deformation, F[0,1]=" + str(F[0, 1]))
	print("TesselationWrapper OK")
