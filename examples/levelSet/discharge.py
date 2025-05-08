# -*- encoding=utf-8 -*-
# jerome.duriez@inrae.fr
# Simulation script for the discharge simulation discussed in the Section 6 of [Duriez2021b]

from yade import ymport, plot, timing
import time

execfile('ramFP.py')
execfile('spaceRot.py')

O.timingEnabled = True
O.trackEnergy = True

#########################
### Engine definition ###
#########################

O.dt = 25 * 1.e-6
nIter = 56000  # 56000 * 25us = 1.4 s
print('Using O.dt = ', O.dt, '; we should run', nIter, 'iterations')
dataPeriod = 50
O.materials.append(FrictMat(density=2650, frictionAngle=radians(25), label='frictional'))  # kn and kt will be directly assigned below, 2650 silica's density
O.engines = [
        ForceResetter(),
        InsertionSortCollider(
                [Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0
        )  # no need to wait for IScollider to decide himself to use 0. Note that I can not use anything else for now: the bounds have to match the grids, and those do not extend much more than the bodies...
        ,
        InteractionLoop(
                [Ig2_LevelSet_LevelSet_ScGeom(), Ig2_Wall_LevelSet_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=1.e5), ks=MatchMaker(algo='val', val=7.e4))
                ],  # 1e5 N/m necessary to keep overlap minimal ?
                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False, label='csLaw2')],
                label='ILoop'
        ),
        NewtonIntegrator(kinSplit=True, damping=0.3, label='newton', gravity=(0, 0, -9.8)),
        PyRunner(command='saveThings(epP0,totP)', iterPeriod=dataPeriod, label='saveData', initRun=True)
]


def saveThings(epP0, totP):
	if 'plastDissip' in O.energy.keys():
		plastW = O.energy['plastDissip']
	else:
		plastW = 0
	if 'elastPotential' in O.energy.keys():
		elP = O.energy['elastPotential']
	else:
		elP = 0
	if 'nonviscDamp' in O.energy.keys():  # using just **O.energy in plot.addData() would be nice as well...
		dampW = O.energy['nonviscDamp']
	else:
		dampW = 0
	maxZb = idBmax = -1
	for b in O.bodies:
		if b.state.pos[2] > maxZb:
			maxZb = b.state.pos[2]
			idBmax = b.id
	epP = epP0 + O.energy['gravWork']  # weird sign convention for gravWork => "+O.energy['gravWork']"
	sumE = elP + epP + plastW + O.energy['kinTrans'] + O.energy['kinRot'] + dampW
	fOnWall = O.forces.f(floor)[2]  # should equal totP once everyone is at equilibrium
	fOnWallRel = fOnWall / totP
	plot.addData(
	        it=O.iter,
	        fOnWallRel=fOnWallRel,
	        elE=elP,
	        kinT=O.energy['kinTrans'],
	        kinR=O.energy['kinRot'],
	        dampW=dampW,
	        epP=epP,
	        fricW=plastW,
	        maxZb=maxZb,
	        idBmax=idBmax,
	        sumE=sumE,
	        speed=O.speed,
	        tILoop=ILoop.execTime * 1.e-3  # now in us
	        ,
	        time=O.time,
	        uF=unbalancedForce(),
	        zc=avgNumInteractions()
	)


plot.plots = {'it': ('kinT', 'elE', 'fricW', 'dampW', 'epP', 'sumE'), 'it ': ('zc', 'fOnWallRel'), ' it ': 'uF', ' it': ('kinT', 'kinR')}

############################
### Particles definition ###
############################

nB = 1000
rcar = 0.01  # longest half-length of se
epsVal = numpy.array([[0.1, 0.5], [0.1, 1], [1, 0.5], [1.4, 1.2], [0.4, 1.6]])
rVal = numpy.array([[0.58, 1, 0.83], [0.42, 1, 0.83], [0.42, 1, 0.83], [0.5, 0.7, 1.], [0.4, 1., 0.8]])
nN, prec = 2000, 20  # LS-DEM choices, where nN = number of boundary nodes and prec is the grid resolution with respect to the min(length)

##################################
### Initial packing parameters ###
##################################


def heightCloud(base, nbr, rad, poros):
	'''Height of a cloud volume of a given (square) *base* length, enclosing *nbr* spheres of radius *rad*, with a *poros* porosity'''
	return (nbr * 4 / 3. * pi * rad**3 / ((1. - poros) * base**2))


rSphCloud = rcar * max(
        numpy.linalg.norm(rVal[o]) for o in range(rVal.shape[0])
)  # particles of the cloud radii have to be sqrt(sum(radii)^2) apart to be sure they (their circumscribed circle) do not touch the neighbour. And we take the maximum of that sum/norm
baseBox = 0.25
baseCloud = baseBox / 1.1
cloudPoros = 0.74  # expected obtained porosity after makeCloud()
h0 = rSphCloud
hCloud = heightCloud(baseCloud, nB, rSphCloud, cloudPoros)
vCloud = hCloud * baseCloud**2

sp = SpherePack()
sphFile = 'discharge.spheres'
if os.path.exists(sphFile):
	sp.load(sphFile)
	print('Existing cloud just loaded')
else:
	sp.makeCloud(rMean=rSphCloud, minCorner=(-baseCloud / 2., -baseCloud / 2., h0), maxCorner=(baseCloud / 2., baseCloud / 2., h0 + hCloud), num=nB)
	while len(sp) != nB:
		print('Giving it another try to obtain', nB, 'bodies')
		sp = SpherePack()
		sp.makeCloud(rMean=rSphCloud, minCorner=(-baseCloud / 2., -baseCloud / 2., h0), maxCorner=(baseCloud / 2., baseCloud / 2., h0 + hCloud), num=nB)
	print('New cloud just build')
	sp.save(sphFile)

########################################
### Creating the DE as per the above ###
########################################

nGP, nNtot = 0, 0  # nGP = nNtot = 0 may as well work for an initialization wo binding together destinies for scalars
epP0, totP, volS = 0, 0, 0
grids, distFields = [None] * 5, [None] * 5  # with grids = distFields = .. grids will always be the same than distFields..
ramFPini = ramFP()
tStart = time.time()
for shape in range(5):  # establishing the list of 5 possible distField once for all
	dimsP = rVal[shape]
	rx, ry, rz = dimsP[0] * rcar, dimsP[1] * rcar, dimsP[2] * rcar
	epsE, epsN = epsVal[shape][0], epsVal[shape][1]
	lsSe = levelSetBody("superellipsoid", extents=(rx, ry, rz), epsilons=(epsE, epsN), spacing=2 * min(rx, ry, rz) / prec, nSurfNodes=2, nodesPath=2)
	grids[shape] = lsSe.shape.lsGrid
	distFields[shape] = lsSe.shape.distField
print('5 (grid ; distance field) pairs computed in', time.time() - tStart, 's and', ramFP() - ramFPini, 'MB of RAM')
for sph in sp:  # before directly using for the present definition of bodies
	currB = len(O.bodies)
	dimsP = rVal[currB % 5]
	O.bodies.append(
	        levelSetBody(
	                center=(sph[0][0], sph[0][1], sph[0][2]),
	                grid=grids[currB % 5],
	                distField=distFields[currB % 5],
	                nSurfNodes=nN,
	                nodesPath=2,
	                dynamic=True,
	                orientation=spaceRot(currB, nB)
	        )
	)
	b = O.bodies[currB]
	nGP += b.shape.lsGrid.nGP.prod()
	nNtot += len(b.shape.surfNodes)
	epP0 += b.state.mass * newton.gravity.norm() * b.state.pos[2]
	totP += b.state.mass * newton.gravity[2]
	volS += b.shape.volume()
tEnd = time.time()
print(
        len(O.bodies), 'LS superellipsoids for a total of', nGP, 'grid points and', nNtot, 'boundary nodes generated in', tEnd - tStart, 's, and',
        ramFP() - ramFPini, 'MB of RAM'
)
print('And a solid volume of', volS, 'i.e. a se cloud porosity of', (vCloud - volS) / vCloud)

floor = O.bodies.append(wall((0, 0, 0), 2))  # ground wall will have the same frictional material, and be non dynamic
O.materials.append(FrictMat(frictionAngle=0, label='frictionless'))
O.bodies.append(wall((-baseBox / 2, 0, 0), 0))
O.bodies.append(wall((baseBox / 2, 0, 0), 0))
O.bodies.append(wall((0., -baseBox / 2, 0), 1))
O.bodies.append(wall((0., baseBox / 2, 0), 1))

print('Initial state has', epP0, 'gravitational energy')

###########
### Run ###
###########

O.run()  # or O.run(nIter+1,True)
plot.plot()
