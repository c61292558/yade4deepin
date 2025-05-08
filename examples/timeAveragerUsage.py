# 2024 © Benjamin Dedieu <benjamin.dedieu@proton.me>
"""
This example illustrates the use of class `TimeAverager` to retrieve averaged data 
of position, velocity force and contact force repartition for 2 particles called the 
intruders entrained in a gravity-driven dry mono-disperse granular flow on an inclined plane.
"""

from yade import pack, qt
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

############
# PARAMETERS

diameterPart = 2e-3  # Diameter of the bed particles [m]
diameterPartI = 4e-3  # Diameter of the intruders [m]
nbIntruders = 2  # number of intruders
positionI = 4 * diameterPart  # target intruder center position from the ground
densPart = 2500  # Density of particles [kg/m3]
phiMax = 0.61  # Max value of the dense packing solid volume fraction
phiMean = 0.6  # Solid volume fraction only used to compute the number of beads
# NbBeads = phiMean*height*width*length/Vpart
length = 15 * diameterPart  # Length of the domain [m]
width = 8 * diameterPart  # Width of the domain [m]
height = 6 * diameterPart  # Height of the domain [m]
slopeAngle = pi / 6  # Slope angle of the inclined plane [rad]
gravityVector = Vector3(9.81 * sin(slopeAngle), 0.0, -9.81 * cos(slopeAngle))
dtSave = 0.1  # time step for saving data [s]
simuTime = 5  # simulation time after gravity deposition [s]
gridDensity = 12  # angle between two points of the grid to compute the
# contact force field on the intruder surface, in degrees
gaussianKernelWidth = 2  # gaussian parameter to diffuse the contact force
# on the grid, normalized by the distance betweeen
# 2 points on the grid

# Estimated max particle pressure from the static load
maxPressure = densPart * phiMax * height * 9.81

# Evaluate the minimal normal stiffness to be in the rigid particle limit
# (cf Roux and Combe 2002)
normalStiffness = maxPressure * diameterPart * 1e4

# Young modulus of the particles from the stiffness wanted.
young = normalStiffness / diameterPart

# Create materials
O.materials.append(ViscElMat(en=0.5, et=0., young=young, poisson=0.5, density=densPart, frictionAngle=atan(0.4), label='Mat'))

# Define grid points on the surface of the intruder to compute stress field from
# contact points. The grid is in cartesian coordinates (simpler than spherical
# coordinates for computing distance to contact positions)
# The grid is generated using the Fibonacci lattice algorythm
Ngrid = int(np.ceil(4 * np.pi / (np.deg2rad(gridDensity)**2)))
grid = np.empty((Ngrid, 3))
goldenRatio = (1 + 5**0.5) / 2
for i in range(Ngrid):
	# compute theta and phi, the sperical coordinates
	theta = 2 * np.pi * i / goldenRatio
	phi = np.arccos(1 - 2 * (i + 0.5) / Ngrid)  # the `+0.5` avoids to have singularity
	# points at the poles
	# convert to cartesian coordinates
	grid[i, :] = diameterPartI / 2 * np.array([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
# compute sigma for the gaussian diffusion of contact forces:
sigma = gaussianKernelWidth * diameterPartI * gridDensity * pi / 180

####################
# FRAMEWORK CREATION

# Compute number of particles to fit in the height
numberPart = int(phiMean * height * length * width / pi / diameterPart**3 * 6)

# Define deposition height for the particle cloud
depositionHeight = 4 * height

# Define number of particles and deposition height for cloud of particles below and above the intruders
heightBelowI = positionI - 0.5 * diameterPartI
heightAboveI = height - heightBelowI
depositionHeightBelowI = heightBelowI / height * depositionHeight
depositionHeightAboveI = depositionHeight - depositionHeightBelowI
numberPartBelowI = int(depositionHeightBelowI / depositionHeight * numberPart)
numberPartAboveI = numberPart - numberPartBelowI
depositionHeightI = 1.2 * diameterPartI
depositionHeightTot = depositionHeightBelowI + depositionHeightI + depositionHeightAboveI

# Define reference positions of the inclined plane
bottomPos = 2 * diameterPart

# Create a periodic box (larger than the flume to avoid bugs):
O.periodic = True
boxHeight = bottomPos + depositionHeightTot
O.cell.setBox(length, width, boxHeight)

# Create walls
bottomWall = box(
        center=(length / 2.0, width / 2.0, bottomPos),
        extents=(1e6, 1e6, 0),
        fixed=True,
        wire=False,
        color=(0, 1, 0),
        material='Mat',
)
O.bodies.append(bottomWall)

# Create the particle cloud
cloudBelowI = pack.SpherePack()
cloudBelowI.makeCloud(
        minCorner=(0, 0, bottomPos + 1e-4),
        maxCorner=(length, width, bottomPos + depositionHeightBelowI),
        rMean=diameterPart / 2,
        num=numberPartBelowI,
)
cloudBelowI.toSimulation(material='Mat', color=(0, 0, 1))
cloudAboveI = pack.SpherePack()
cloudAboveI.makeCloud(
        minCorner=(0, 0, bottomPos + depositionHeightBelowI + depositionHeightI),
        maxCorner=(length, width, depositionHeightTot),
        rMean=diameterPart / 2,
        num=numberPartAboveI,
)
cloudAboveI.toSimulation(material='Mat', color=(0, 0, 1))

# Define left and right Limit of intruder zone :
leftLimitI = diameterPartI / 2
rightLimitI = diameterPartI / 2

# Create the intruders (equaly distributed along x axis)
intruders = []
for i in range(nbIntruders):
	O.bodies.append(
	        sphere(
	                center=((i + 0.5) * length / nbIntruders, leftLimitI, bottomPos + depositionHeightBelowI + depositionHeightI / 2),
	                radius=diameterPartI / 2.,
	                color=(1, 0, 0),
	                material='Mat',
	                wire=False
	        )
	)
	intruders.append(O.bodies[-1])

# Evaluate the deposition time considering the free-fall time of the highest particle to the ground
depositionTime = sqrt(depositionHeightTot * 2 / abs(gravityVector[2]))

###################################
# ENGINES AND PY-RUNNER DEFINITION

O.engines = [
        ForceResetter(),
        InsertionSortCollider(
                [Bo1_Sphere_Aabb(), Bo1_Box_Aabb()],
                allowBiggerThanPeriod=True,
        ),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
                [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
                [Law2_ScGeom_ViscElPhys_Basic()],
        ),
        PyRunner(
                command='initialize()',
                virtPeriod=0.01,
                label='init',
        ),
        TimeAverager(
                ids=[intr.id for intr in intruders],
                computeContactForceField=True,
                sigma=sigma,
                grid=grid,
                label='timeAverager',
                dead=True,
        ),
        PyRunner(
                command='saveData()',
                virtPeriod=dtSave,
                label='dataSaver',
                dead=True,
        ),
        GlobalStiffnessTimeStepper(
                defaultDt=1e-4,
                viscEl=True,
                timestepSafetyCoefficient=0.7,
        ),
        NewtonIntegrator(gravity=gravityVector, damping=0.2),
]


def initialize():
	"""Wait for gravity deposition then start TimeAverager and saveData and 
    initialize quantities"""
	global step, time, Nstep
	if O.time >= depositionTime:
		init.dead = True
		timeAverager.dead = False
		dataSaver.dead = False
		Nstep = int(simuTime / dtSave)
		for intr in intruders:
			intr.instantPos = np.empty((Nstep, 3))
			intr.averagedPos = np.empty((Nstep, 3))
			intr.instantVel = np.empty((Nstep, 3))
			intr.averagedVel = np.empty((Nstep, 3))
			intr.instantForce = np.empty((Nstep, 3))
			intr.averagedForce = np.empty((Nstep, 3))
			intr.averagedContactForceField = np.empty((Nstep, Ngrid, 3))
		time = np.empty((Nstep))
		step = 0
		timeAverager.initialization()


def saveData():
	global step, time, cellArea, FcFieldMean, pressureFieldMean, fig, fig0, fig1, fig2
	if step <= Nstep - 1:
		time[step] = O.time
		for k, intr in enumerate(intruders):
			# Get instantaneous quantities direclty from body objects
			intr.instantPos[step] = intr.state.pos
			intr.instantVel[step] = intr.state.vel
			intr.instantForce[step] = O.forces.f(intr.id)
			# Get averaged quantities with timeAverager.get... methods
			intr.averagedPos[step] = timeAverager.getPos(intr.id)
			intr.averagedVel[step] = timeAverager.getVel(intr.id)
			intr.averagedForce[step] = timeAverager.getForce(intr.id)
			intr.averagedContactForceField[step] = timeAverager.getContactForceField(intr.id)

		# update step and re-initialize timeAverager
		print(f"Data saved for step {step}/{Nstep-1}")
		timeAverager.initialization()
		step += 1

	if step == Nstep:
		# Pause the simulation
		O.pause()

		# Prepare plots
		plt.ion()
		cmap = plt.get_cmap("viridis", nbIntruders)

		# Plot intruder z position over time for each intruder
		fig0, ax0 = plt.subplots(constrained_layout=True)
		for k, intr in enumerate(intruders):
			ax0.plot(time, intr.averagedPos[:, 2] * 1000, '-', color=cmap(k), label="intruder " + str(k) + " (averaged)")
			ax0.plot(time, intr.instantPos[:, 2] * 1000, '-', color=cmap(k), alpha=0.2, label="intruder " + str(k) + " (intantaneous)")
		ax0.set_xlabel("$time$ [s]")
		ax0.set_ylabel("$z$ [mm]")
		ax0.legend()

		# Plot intruder x velocity over time for each intruder
		fig1, ax1 = plt.subplots(constrained_layout=True)
		for k, intr in enumerate(intruders):
			ax1.plot(time, intr.averagedVel[:, 0], '-', color=cmap(k), label="intruder " + str(k) + " (averaged)")
			ax1.plot(time, intr.instantVel[:, 0], '-', color=cmap(k), alpha=0.2, label="intruder " + str(k) + " (intantaneous)")
		ax1.set_xlabel("$time$ [s]")
		ax1.set_ylabel("$v_x$ [m/s]")
		ax1.legend()

		# Plot intruder z total force over time for each intruder
		fig2, ax2 = plt.subplots(constrained_layout=True)
		for k, intr in enumerate(intruders):
			ax2.plot(time, intr.averagedForce[:, 2], '-', color=cmap(k), label="intruder " + str(k) + " (averaged)")
			ax2.plot(time, intr.instantForce[:, 2], '-', color=cmap(k), alpha=0.2, label="intruder " + str(k) + " (intantaneous)")
		ax2.set_xlabel("$time$ [s]")
		ax2.set_ylabel("$F_z$ [N]")
		ax2.legend()

		# Plot pressure repartition at the surface of the intruder 0
		FcFieldMean = np.mean(intruders[0].averagedContactForceField, axis=0)
		# compute the normal force
		FcnFieldMean = np.empty((Ngrid))
		for i in range(Ngrid):
			FcnFieldMean[i] = np.dot(grid[i] / diameterPartI, -FcFieldMean[i])
		# Adimension by the max of normal force
		FcnFieldMean /= max(FcnFieldMean)
		fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, tight_layout=True)
		x, y, z = grid[:, 0] * 1000, grid[:, 1] * 1000, grid[:, 2] * 1000
		norm = mpl.colors.Normalize(vmin=min(FcnFieldMean), vmax=max(FcnFieldMean))
		cmap = plt.get_cmap("viridis")
		colors = np.array([cmap(norm(F)) for F in FcnFieldMean])
		ax.scatter(x, y, z, color=colors)
		ax.set_xlabel("x [mm]")
		ax.set_ylabel("y [mm]")
		ax.set_zlabel("z [mm]")
		fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, aspect=60, pad=0.1, label="P / Pmax [-]")
		fig.suptitle("Pressure repartion on Intruder 0")


###################
# RUN WITH QT VIEW
view = qt.View()
view.viewDir = (0, 1, 0)
view.fitAABB([0, 0, 0], [length, width, boxHeight])
O.run()
