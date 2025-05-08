if ('LS_DEM' in features):
	# Handy functions for comparing obtained quantities with expected ones (TODO: this is a copy-paste from checkLS_DEM.py...)
	######################################################################
	def equalVectors(n, nRef, tol=0.02):
		'''Whether n and nRef are equal or opposite vectors, up to a relative tolerance tol (that applies to the norm)'''
		if (n - nRef).norm() > tol * nRef.norm() and (n + nRef).norm() > tol * nRef.norm():
			return False
		else:
			return True

	def equalNbr(x, xRef, tol=0.02):
		'''Whether x and xRef are equal numbers, up to a relative tolerance tol'''
		if abs(x - xRef) / abs(xRef) > tol:
			return False
		else:
			return True

	# Let us check MultiScGeom behavior for a LStwin of a clump on a ground
	#######################################################################
	# NB: because of the LStwin shape, physics is actually in a bifurcation regime and the exact dynamics may vary (e.g., on compilation options)

	for repeat in range(2):  # we will repeat twice ~ the same thing, with either a LSbox or a true Wall for the ground
		# The bodies
		############
		O.bodies.appendClumped([sphere((0, 0, z), 1) for z in [-0.5, 0.5]])  # a 3*1*1 long body
		clumpB = O.bodies[2]
		lsClumpId = O.bodies.append(
		        levelSetBody(clump=clumpB.shape, center=(1.6, 0, 0), nSurfNodes=3800, spacing=0.1)
		)  # highest dimension along x actually
		lsClump = O.bodies[lsClumpId]
		O.bodies.erase(clumpB.id, True)

		#now the two versions of the x=0 ground:
		if repeat == 0:  # using the LSbox (jduriez note: this is aor8)
			hBox = 0.5
			posB = Vector3(-hBox, 0, 0)
			bId = O.bodies.append(levelSetBody("box", extents=(hBox, 5, 5), spacing=0.25, center=posB, dynamic=False))
		else:  # using a true Wall (jduriez note: this is aor8wall)
			bId = O.bodies.append(wall(0, 0))

		# The engines
		#############
		knVal = 2.4e8
		ksVal = 0.5 * knVal
		O.engines = [
		        ForceResetter(),
		        InsertionSortCollider([Bo1_LevelSet_Aabb()] if repeat == 0 else [Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0),
		        InteractionLoop(
		                [Ig2_LevelSet_LevelSet_MultiScGeom()] if repeat == 0 else [Ig2_Wall_LevelSet_MultiScGeom()],
		                [Ip2_FrictMat_FrictMat_MultiFrictPhys(kn=knVal, ks=ksVal)],
		                [Law2_MultiScGeom_MultiFrictPhys_CundallStrack(sphericalBodies=False)]
		        ),
		        NewtonIntegrator(gravity=Vector3(-9.8, 0, 0), label='ni', damping=0.3)
		]
		O.dt = 0.7 * (
		        lsClump.shape.volume() * lsClump.mat.density / knVal
		)**0.5  # 0.8 * .. ~ the maximum value for some numeric settings and above maximum for others. We use 0.7
		O.run(2500, True)  # 2200 ~ the smallest value with O.dt = 0.8 * .. possible to lead to equilibrium. Using 2500 with 0.7 * ..

		if not O.interactions.has(bId, lsClump.id, True):
			raise YadeCheckError(
			        "We're missing the", O.bodies[bId].shape, "-", lsClump.shape, "expected interaction (while being in repeat =", repeat,
			        "step of the script)"
			)
		cont = O.interactions[bId, lsClump.id]
		weight = lsClump.shape.volume() * lsClump.mat.density * ni.gravity
		if not equalVectors(O.forces.f(bId), weight, 0.05):
			raise YadeCheckError("Incorrect equilibrium state between applied force onto the ground =", O.forces.f(bId), "and weight", weight)
		expectedCtctPts = [
		        2, 3, 4
		]  # while 2 is the true answer, imperfect script (and variable compilation) conditions may induce 1 or 2 extra contacting nodes
		obtainedCtctPts = len(cont.geom.contacts)
		if not obtainedCtctPts in expectedCtctPts:
			raise YadeCheckError(
			        "Changed behavior in MultiScGeom, got :", len(cont.geom.contacts), "vs", str(expectedCtctPts[0]),
			        "expected (with a +2 tolerance)"
			)
		O.reset()
else:
	print("Skip checkMultiScGeom, LS-DEM feature not available")
