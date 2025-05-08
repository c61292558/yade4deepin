################################################################################
# # Checking viscoelastic boundary condition: VESupportEngine()
# A sphere is loaded with constant a force. Measured displacement is compared to the result of the analytical solution.
# Initially there is a difference between plots because the analytical solution is simplified (ommits mass).
# Hence, only the displacements measured and computed at the end of the simulations are compared.
################################################################################

# PARAMETERS
k1 = 500
k2 = 100
c1 = 1000
c2 = 50
force = 100
# BODIES

sph = O.bodies.append(sphere(center=(1, 0, 0), radius=0.1))
O.forces.setPermF(sph, (0, 0, force))

# ENGINES

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        VESupportEngine(bIds=[sph], k1=k1, k2=k2, c1=c1, c2=c2),
        NewtonIntegrator(damping=0.0, label='newton'),
]

########
O.dt = 0.8 * PWaveTimeStep()

N = int(5 / O.dt)

O.run(N, wait=True)

# Comparison (if ok, YadeCheckError Exception is not thrown)
tolerance = 0.001

d1 = O.bodies[sph].state.pos[2]
# The reference value of d2 can be computed as follows:
d2 = force / k1 + force * O.time / c1 + (force / k2) * (1 - yade.math.exp(-k2 * O.time / c2))

if (abs(d1 - d2 > tolerance)):
	raise YadeCheckError("abs(d1-d2 > tolerance")
