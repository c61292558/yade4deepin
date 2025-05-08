"""
This example illustrates the use of the function `getSlicedProfiles` to get average
profiles of particle volume fraction and particle velocity over a given time and along
a specified axis.

The configuration of this example is a gravity-driven dry mono-disperse granular
flow in a flume (inclined plane with side walls).

Three examples of profiles are given :
- example 1 : vertical (z) profiles at the center of the flume
- example 2 : vertical (z) profiles at the side walls
- example 3 : cross (y) profiles between 2 elevations
"""

from yade import pack, qt
import numpy as np
from matplotlib import pyplot as plt

############
# PARAMETERS

diameterPart = 2e-3  # Diameter of the particles [m]
densPart = 2500  # Density of particles [kg/m3]
phiPartMax = 0.6  # Max particle concentration estimation
phiPartMean = 0.58  # Mean particle concentration estimation
length = 15 * diameterPart  # Length of the flume [m]
width = 8 * diameterPart  # Width of the flume [m]
height = 6 * diameterPart  # Height of the flume [m]
slopeAngle = pi / 4  # Slope angle of the flume [rad]
delta = diameterPart / 30  # Discretisation step for the profiles
gravityVector = Vector3(9.81 * sin(slopeAngle), 0.0, -9.81 * cos(slopeAngle))
dtSave = 0.1  # time step for saving profiles [s]
timeAverage = 3  # time during which we average profiles after deposition [s]

# Define slice elevation limits for example 3
zBottom = 3 * diameterPart
zTop = 5 * diameterPart

# Estimated max particle pressure from the static load
maxPressure = densPart * phiPartMax * height * 9.81

# Evaluate the minimal normal stiffness to be in the rigid particle limit
# (cf Roux and Combe 2002)
normalStiffness = maxPressure * diameterPart * 1e4

# Young modulus of the particles from the stiffness wanted.
young = normalStiffness / diameterPart

# Create materials
O.materials.append(ViscElMat(en=0.5, et=0., young=young, poisson=0.5, density=densPart, frictionAngle=atan(0.4), label='Mat'))

####################
# FRAMEWORK CREATION

# Compute number of particles to fit in the height
numberPart = int(phiPartMean * height * length * width / pi / diameterPart**3 * 6)

# Define deposition height for the particle cloud
depositionHeight = 4 * height

# Define reference positions of the flume
bottomPos = 2 * diameterPart
leftPos = diameterPart

# Create a periodic box (larger than the flume to avoid bugs):
O.periodic = True
boxWidth = width + 2 * leftPos
boxHeight = bottomPos + depositionHeight
O.cell.setBox(length, boxWidth, boxHeight)

# Create walls
bottomWall = box(
        center=(length / 2.0, boxWidth / 2.0, bottomPos),
        extents=(1e6, boxWidth / 2.0, 0),
        fixed=True,
        wire=False,
        color=(0, 1, 0),
        material='Mat',
)
leftWall = box(
        center=(length / 2.0, leftPos, boxHeight / 2.0),
        extents=(1e6, 0, boxHeight / 2.0),
        fixed=True,
        wire=True,
        color=(1, 0, 0),
        material='Mat',
)
rightWall = box(
        center=(length / 2.0, leftPos + width, boxHeight / 2.0),
        extents=(1e6, 0, boxHeight / 2.0),
        fixed=True,
        wire=True,
        color=(1, 0, 0),
        material='Mat',
)
O.bodies.append([bottomWall, leftWall, rightWall])

# Create the particle cloud
epsilon = diameterPart / 4  # to avoid bugs with initial contacts between particles and walls
cloud = pack.SpherePack()
cloud.makeCloud(
        minCorner=(0, leftPos + epsilon, bottomPos + 1e-4),
        maxCorner=(length, leftPos + width - epsilon, boxHeight),
        rMean=diameterPart / 2,
        num=numberPart,
)
cloud.toSimulation(material='Mat', color=(0, 0, 1))

# Evaluate the deposition time considering the free-fall time of the highest particle to the ground
depositionTime = sqrt(depositionHeight * 2 / abs(gravityVector[2]))

#################
# INITIALIZATION

# Compute y and z axis for latter plots
yAxis = np.array([i * delta / diameterPart for i in range(int(width / delta))])
zAxis = np.array([i * delta / diameterPart for i in range(int(height / delta))])

# Initialize arrays of concentration and velocity for the 3 examples
Nstep = int(timeAverage / dtSave)
phiPart1 = np.empty([Nstep, len(zAxis)])
vPart1 = np.empty([Nstep, len(zAxis), 3])
phiPart2 = np.empty([Nstep, len(zAxis)])
vPart2 = np.empty([Nstep, len(zAxis), 3])
phiPart3 = np.empty([Nstep, len(yAxis)])
vPart3 = np.empty([Nstep, len(yAxis), 3])

# Initialize the step to fill the phiPart and vPart array
step = 0

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
        GlobalStiffnessTimeStepper(
                defaultDt=1e-4,
                viscEl=True,
                timestepSafetyCoefficient=0.7,
        ),
        NewtonIntegrator(gravity=gravityVector, damping=0.2),
        PyRunner(command='saveProfiles()', virtPeriod=dtSave),
]


def saveProfiles():
	global phiPart1, vPart1, phiPart2, vPart2, phiPart3, vPart3, step

	if O.time >= depositionTime and O.time < depositionTime + timeAverage:

		# example 1 : get z profiles at the center of the flume on a slice large
		# of `2*diameterPart` in y-direction, with `int(height/delta)` intervals
		# of discretisation of size `delta` in z-direction.
		phiPart1[step, :], vPart1[step, :, :] = getSlicedProfiles(
		        vCell=2 * diameterPart * length * delta,
		        nCell=int(height / delta),
		        dP=delta,
		        refP=bottomPos,
		        refS=leftPos,
		        dirS=1,
		        dirP=2,
		        sliceCenters=[width / 2],
		        sliceWidths=[2 * diameterPart]
		)

		# example 2 : get z profiles at the box walls on two sub-slices large of
		# `diameterPart` in y-direction, with `int(height/delta)` intervals of
		# discretisation of size `delta` in z-direction.
		phiPart2[step, :], vPart2[step, :, :] = getSlicedProfiles(
		        vCell=2 * diameterPart * length * delta,
		        nCell=int(height / delta),
		        dP=delta,
		        refP=bottomPos,
		        refS=leftPos,
		        dirS=1,
		        dirP=2,
		        sliceCenters=[diameterPart / 2, width - diameterPart / 2],
		        sliceWidths=[diameterPart] * 2
		)

		# example 3 : get y profiles on a slice which z position is between 2 and
		# 4 diameters from the bottom of the flume with `int(width/delta)`
		# intervals of discretisation of size `delta` in y-direction.
		phiPart3[step, :], vPart3[step, :, :] = getSlicedProfiles(
		        vCell=(zTop - zBottom) * length * delta,
		        nCell=int(width / delta),
		        dP=delta,
		        refP=leftPos,
		        refS=bottomPos,
		        dirS=2,
		        dirP=1,
		        sliceCenters=[(zTop + zBottom) / 2],
		        sliceWidths=[zTop - zBottom]
		)

		# update step
		print(f"Profiles saved for step {step}/{Nstep}.")
		step += 1

	if O.time >= depositionTime + timeAverage:
		# Pause the simulation
		O.pause()

		# Average Profiles over all the time steps
		phiPart1 = np.mean(phiPart1, axis=0)
		vPart1 = np.mean(vPart1, axis=0)
		phiPart2 = np.mean(phiPart2, axis=0)
		vPart2 = np.mean(vPart2, axis=0)
		phiPart3 = np.mean(phiPart3, axis=0)
		vPart3 = np.mean(vPart3, axis=0)

		# Extract streamwise velocities
		vxPart1 = vPart1[:, 0]
		vxPart2 = vPart2[:, 0]
		vxPart3 = vPart3[:, 0]

		# Plot z profiles (examples 1 and 2)
		fig, axs = plt.subplots(1, 2, sharey=True)
		fig.suptitle("z profiles (Examples 1 and 2)")
		axs[0].plot(phiPart1, zAxis, label="Example 1 : at the center")
		axs[0].plot(phiPart2, zAxis, label="Example 2 : at the walls")
		axs[1].plot(vxPart1, zAxis, label="Example 1 : at the center")
		axs[1].plot(vxPart2, zAxis, label="Example 2 : at the walls")
		axs[0].set_xlabel(r"$\phi$ [-]")
		axs[1].set_xlabel("$v_x$ [m/s]")
		axs[0].set_ylabel("$z$ [diameters]")
		axs[1].legend()
		fig.show()

		# Plot y profiles (example 3)
		fig, axs = plt.subplots(1, 2, sharey=True)
		fig.suptitle("y cross profiles (Example 3)")
		axs[0].plot(phiPart3, yAxis)
		axs[1].plot(vxPart3, yAxis)
		axs[0].set_xlabel(r"$\phi$ [-]")
		axs[1].set_xlabel("$v_x$ [m/s]")
		axs[0].set_ylabel("$y$ [diameters]")
		fig.show()


###################
# RUN WITH QT VIEW
view = qt.View()
view.viewDir = (0, 1, 0)
view.fitAABB([0, 0, 0], [length, boxWidth, boxHeight])
O.run()
