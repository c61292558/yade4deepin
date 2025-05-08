import random
import math as m

random.seed(11)

# Create a sphere cloud
sp = pack.SpherePack()
sp.makeClumpCloud(
        (0, 0, 0), (10, 10, 10),
        [pack.SpherePack([((-0.04, -0.04, -0.04), 0.02 * sqrt(3)), ((0.0, 0.0, 0.0), 0.02 * sqrt(3)), ((0.04, 0.04, 0.04), 0.02 * sqrt(3))])],
        periodic=False,
        num=100,
        seed=11
)
sp.toSimulation()

# Ensure that all particles have the same kinetic energy
Vel = 2.0
for b in O.bodies:
	if isinstance(b.shape, Clump):
		theta = random.uniform(0, m.pi)
		phi = random.uniform(0, m.pi)
		vv = Vector3(m.sin(theta) * m.cos(phi), m.sin(theta) * m.sin(phi), m.cos(theta))
		R = b.state.ori.toRotationMatrix()
		b.state.vel = 0.05 * Vel * vv
		b.state.angVel = Vel * vv
		b.state.angMom = (R * b.state.inertia.asDiagonal() * R.transpose()) * b.state.angVel

# Set time step
O.stopAtTime = 2.0
O.dt = 0.0001

# Engines
O.trackEnergy = True
O.engines = [
        ForceResetter(),
        InsertionSortCollider([]),
        InteractionLoop([], [], []),
        NewtonIntegrator(damping=0.0, gravity=[0, 0, 0], label='newton', kinSplit=True, exactAsphericalRot=True)
]

O.step()
E0_1 = utils.kineticEnergy()
E0_2 = O.energy['kinRot'] + O.energy['kinTrans']

if abs(E0_1 - E0_2) / E0_1 > 1.5e-2:
	raise YadeCheckError("1. Energy calculated in utils.kineticEnergy (", E0_1, ") does not match with O.trackEnergy (", E0_2, ").")

O.run()
O.wait()

Ef_1 = utils.kineticEnergy()
Ef_2 = O.energy['kinRot'] + O.energy['kinTrans']

if abs(Ef_1 - Ef_2) / Ef_1 > 1.5e-2:
	raise YadeCheckError("2. Energy calculated in utils.kineticEnergy (", Ef_1, ") does not match with O.trackEnergy (", Ef_2, ").")

err1 = abs(E0_1 - Ef_1) / E0_1
err2 = abs(E0_2 - Ef_2) / E0_2

if err1 > 3e-2 or err2 > 3e-2:
	raise YadeCheckError("Energy relative error over time is too big (", err1, " or ", err2, ")")
