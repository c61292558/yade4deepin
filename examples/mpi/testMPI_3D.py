# -*- encoding=utf-8 -*-

# Possible executions of this script
### Parallel:
# mpiexec -n 4 yade-mpi -n -x testMPIxNxM.py
# mpiexec -n 4 yade-mpi  -n -x testMPIxN.py n # n is the number of workers you wish to apply to the domain.
# yade-mpi -n -x testMPIxN.py
# yade-mpi -n -x testMPIxN.py n
# in last line the optional argument 'n' has the same meaning as with mpiexec
### Openmp:
# yade-mpi -j4 -n -x testMPIxN.py n
### Nexted MPI * OpenMP
# needs testing...
'''
This script simulates spheres falling on a plate using a distributed memory approach based on mpy module. 
The script decomposes a preset 3D domain into n number of subdomains i.e. it is meant to test the strong scalability of yade-mpi
The master process (rank=0) has no spheres assigned; it is in charge of getting the total force on the plate

'''

NSTEPS = 1000  #turn it >0 to see time iterations, else only initilization

#################
# Check MPI world
# This is to know if it was run with or without mpiexec (see preamble of this script)
import os

rank = os.getenv('OMPI_COMM_WORLD_RANK')
if rank is not None:  #mpiexec was used
	rank = int(rank)
	numThreads = int(os.getenv('OMPI_COMM_WORLD_SIZE'))
else:  #non-mpi execution, numThreads will still be used as multiplier for the problem size (2 => multiplier is 1)
	numThreads = 2 if len(sys.argv) < 1 else (int(sys.argv[1]))
	print "numThreads", numThreads

############  Build a scene (we use Yade's pre-filled scene)  ############

from yade import pack
import numpy as np
#mn,mx=Vector3(-10,-10,-10),Vector3(10,10,10) # uniform example
mn, mx = Vector3(-10, -15, -20), Vector3(10, 15, 20)  # heterogeneous geometry example
pred = pack.inAlignedBox(mn, mx)
O.bodies.append(pack.regularHexa(pred, radius=0.8, gap=0))

from yade import domaindecomposition as dd

globalDomain = dd.GlobaldomainCloud(
        mn, mx, numThreads - 1
)  # automatically decompose domain into point cloud based subdomains (random subdomain shapes). Best choice for using random numbers of MPI threads
#globalDomain = dd.Globaldomain(mn,mx,numThreads-1) # automatically decompose the domain in organized grid for all directions for desired number of threads. Can use random numbers of MPI threads, but finicky
#globalDomain = dd.Globaldomain(mn,mx,xDecomp=1,yDecomp=1,zDecomp=numThreads-1) # manually decompose the domain in z direction for given number of threads. Ok choice for using random numbers of MPI threads

# sequential grain colors
import colorsys

colorScale = (Vector3(colorsys.hsv_to_rgb(value * 1.0 / numThreads, 1, 1)) for value in range(0, numThreads))
colors = []
for i in range(numThreads):
	colors.append(next(colorScale))

for b in O.bodies:
	if isinstance(b.shape, Sphere):
		point = dd.Point(b.state.pos[0], b.state.pos[1], b.state.pos[2])
		sd = globalDomain.findSubdomain(point)  # use globalDomain to identify body subdomain
		if rank is not None:
			b.subdomain = sd
		b.shape.color = colors[sd]

walls = aabbWalls([Vector3(mn[0] * 2, mn[1], mn[2]), Vector3(mx[0] * 2, mx[1], mx[2])], thickness=0)
wallIds = O.bodies.append(walls)
WALL_ID = wallIds[2]

collider.verletDist = 0.7
collider.targetInterv = 200
newton.gravity = (0, -10, 0)  #else nothing would move
tsIdx = O.engines.index(timeStepper)  #remove the automatic timestepper. Very important: we don't want subdomains to use many different timesteps...
O.engines = O.engines[0:tsIdx] + O.engines[tsIdx + 1:]
#O.dt=0.00002 #this very small timestep will make it possible to run 2000 iter without merging
O.dt = 0.1 * PWaveTimeStep()  #very important, we don't want subdomains to use many different timesteps...


#########  RUN  ##########
def collectTiming():
	created = os.path.isfile("collect_3D.dat")
	f = open('collect.dat', 'a')
	if not created:
		f.write("numThreads mpi omp Nspheres runtime force\n")
	from yade import timing
	f.write(
	        str(globalDomain.numThreads) + " " + str(os.getenv('OMPI_COMM_WORLD_SIZE')) + " " + os.getenv('OMP_NUM_THREADS') + " " + str(len(O.bodies)) +
	        " " + str(timing.runtime()) + " " + str(O.forces.f(WALL_ID)[1]) + "\n"
	)
	f.close()


if rank is None:  #######  Single-core  ######
	O.timingEnabled = True
	O.run(NSTEPS, True)
	#print "num bodies:",len(O.bodies)
	from yade import timing
	timing.stats()
	collectTiming()
	print "num. bodies:", len([b for b in O.bodies]), len(O.bodies)
	print "Total force on floor=", O.forces.f(WALL_ID)[1]
else:  #######  MPI  ######
	#import yade's mpi module
	from yade import mpy as mp
	# customize
	mp.ACCUMULATE_FORCES = True  #trigger force summation on master's body (here WALL_ID)
	mp.VERBOSE_OUTPUT = False
	mp.ERASE_REMOTE = True  #erase bodies not interacting wit a given subdomain?
	mp.OPTIMIZE_COM = True  #L1-optimization: pass a list of double instead of a list of states
	mp.USE_CPP_MPI = True and mp.OPTIMIZE_COM  #L2-optimization: workaround python by passing a vector<double> at the c++ level

	mp.mpirun(NSTEPS, True)  # passing numThreads
	print "num. bodies:", len([b for b in O.bodies]), len(O.bodies)
	if rank == 0:
		mp.mprint("Total force on floor=" + str(O.forces.f(WALL_ID)[1]))
		collectTiming()
	else:
		mp.mprint("Partial force on floor=" + str(O.forces.f(WALL_ID)[1]))
	mp.mergeScene()
	if rank == 0:
		O.save('mergedScene.yade')
	mp.MPI.Finalize()
#exit()
