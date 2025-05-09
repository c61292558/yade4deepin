# encoding: utf-8

# periodic triaxial test simulation
#
# The initial packing is either
#
# 1. random cloud with uniform distribution, or
# 2. cloud with specified granulometry (radii and percentages), or
# 3. cloud of clumps, i.e. rigid aggregates of several particles
#
# The triaxial consists of 2 stages:
#
# 1. isotropic compaction, until sigmaIso is reached in all directions;
#    this stage is ended by calling compactionFinished()
# 2. constant-strain deformation along the z-axis, while maintaining
#    constant stress (sigmaIso) laterally; this stage is ended by calling
#    triaxFinished()
#
# Controlling of strain and stresses is performed via PeriTriaxController,
# of which parameters determine type of control and also stability
# condition (maxUnbalanced) so that the packing is considered stabilized
# and the stage is done.
#

sigmaIso = -1e5

# download input data from yade-data package
os.system("wget -O capillaryfile.txt -nc https://gitlab.com/yade-dev/yade-data/-/raw/main/capillaryFiles/capillaryfile.txt?ref_type=heads&inline=false")

#import matplotlib
#matplotlib.use('Agg')

# generate loose packing
from yade import pack, qt, plot

O.periodic = True
sp = pack.SpherePack()
if 0:
	## uniform distribution
	sp.makeCloud((0, 0, 0), (2, 2, 2), rMean=.1, rRelFuzz=.3, periodic=True)
else:
	## create packing from clumps
	# configuration of one clump
	c1 = pack.SpherePack([((0, 0, 0), .03333), ((.03, 0, 0), .017), ((0, .03, 0), .017)])
	# make cloud using the configuration c1 (there could c2, c3, ...; selection between them would be random)
	sp.makeClumpCloud((0, 0, 0), (2, 2, 2), [c1], periodic=True, num=500)

# setup periodic boundary, insert the packing
sp.toSimulation()

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=1.2)], [Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay()],
                [Law2_ScGeom_FrictPhys_CundallStrack(neverErase=True)]
        ),
        PeriTriaxController(
                label='triax',
                # specify target values and whether they are strains or stresses
                goal=(sigmaIso, sigmaIso, sigmaIso),
                stressMask=7,
                # type of servo-control
                dynCell=True,
                maxStrainRate=(10, 10, 10),
                # wait until the unbalanced force goes below this value
                maxUnbalanced=.1,
                relStressTol=1e-3,
                # call this function when goal is reached and the packing is stable
                doneHook='compactionFinished()'
        ),
        ######################     INSERT CAPILLARITY ENGINE HERE ############################################
        CapillarityEngine(capillaryPressure=1000, liquidTension=10000, dead=True, label="capillarity"),
        ######################################################################################################
        NewtonIntegrator(damping=.2),
        PyRunner(command='addPlotData()', iterPeriod=100),
]
O.dt = .5 * PWaveTimeStep()


def addPlotData():
	plot.addData(
	        unbalanced=unbalancedForce(),
	        i=O.iter,
	        sxx=triax.stress[0],
	        syy=triax.stress[1],
	        szz=triax.stress[2],
	        exx=triax.strain[0],
	        eyy=triax.strain[1],
	        ezz=triax.strain[2],
	        # save all available energy data
	        Etot=O.energy.total(),
	        **O.energy
	)


# enable energy tracking in the code
O.trackEnergy = True

# define what to plot
plot.plots = {
        'i': ('unbalanced',),
        'i ': ('sxx', 'syy', 'szz'),
        ' i': ('exx', 'eyy', 'ezz'),
        # energy plot
        ' i ': (O.energy.keys, None, 'Etot'),
}
# show the plot
plot.plot()


def compactionFinished():
	# set the current cell configuration to be the reference one
	O.cell.trsf = Matrix3.Identity
	# change control type: keep constant confinement in x,y, 20% compression in z
	triax.goal = (sigmaIso, sigmaIso, -.2)
	triax.stressMask = 3
	# enable capillary forces
	capillarity.dead = False
	# allow faster deformation along x,y to better maintain stresses
	triax.maxStrainRate = (1., 1., .1)
	# next time, call triaxFinished instead of compactionFinished
	triax.doneHook = 'triaxFinished()'
	# do not wait for stabilization before calling triaxFinished
	triax.maxUnbalanced = 10


def triaxFinished():
	print('Finished')
	O.pause()


from yade import timing

O.timingEnabled = True
