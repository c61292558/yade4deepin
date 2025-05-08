# Illustration of MultiScGeom and related classes for non-convex LevelSet shapes
# The script simulates a rounded tetrahedron bouncing along z-axis on a Wall-shaped floor
# By jerome.duriez@inrae.fr and sacha.duverger@inrae.fr

from yade import plot

###################
# Defining bodies #
###################

# Handy function to define the same rounded tetrahedron as considered in https://link.springer.com/article/10.1007/s10035-023-01378-z:
import numpy as np


def get_tet_vertices(sphere_radius=3.101e-3, clump_radius=5e-3, pos=np.array([0, 0, 0]), return_spheres=True, sphere_as_tuple=False):
	"""
	Return a list of four (1,3) numpy array representing the vertices of a tetrahedron.
	If return_sphere, the objects in the list are yade spheres (and not numpy arrays). 
	"""
	tet2sph_dist = (clump_radius - sphere_radius) / 2  # distance between the center of the tetrahedron and the center of the spheres

	e = (clump_radius - sphere_radius) * 4 / np.sqrt(6)  # edge length
	h = np.sqrt(6) / 3 * e  # heigth
	Rcirc = np.sqrt(6) / 4 * e  # radius of the circumsphere
	proj_e2base = np.sqrt(Rcirc**2 - (h - Rcirc)**2)  # projection of an edge to the bottom base

	vs = []
	vs.append(pos + np.array([0, Rcirc, 0]))  # vertex 1 (the one on the top of the pyramid)
	vs.append(pos + np.array([0, Rcirc - h, -proj_e2base]))  # vertex 2 (the one behind the pyramid)
	a = np.sqrt(3 * e**2 / 4) - proj_e2base
	vs.append(pos + np.array([-e / 2, Rcirc - h, a]))  # vertex 3 (the one on the left)
	vs.append(pos + np.array([e / 2, Rcirc - h, a]))  # vertex 4 (the one on the right)

	if not return_spheres:
		return vs
	else:
		if sphere_as_tuple:
			return [(c, sphere_radius) for c in vs]
		else:
			return [sphere(c, sphere_radius) for c in vs]


# Defining the Clump template first:
members = get_tet_vertices()
defaultPos = Vector3(0.000010117244255388888, 1.4167583104224976e-6, -5.669863590932203e-6)  # barycenter of above members with discretization = 20
zC = 0.0038  # will also serve as the initial height of the LS body, has to be more than 3.73 mm to avoid initial contact
wantedPos = Vector3(0, 0, zC)
someOri = Quaternion((1, 0, 0), pi / 2)  # this will not be the actual orientation
for member in members:
	member.state.pos = wantedPos + someOri * (member.state.pos - defaultPos)
O.bodies.appendClumped(members, discretization=20)
clumpB = O.bodies[4]
clump = clumpB.shape

# Now defining the LS body from the above Clump template:
zC = 0.0038
spac = 200  # level set grid spacing, in microns
nN = 1600  # desired number of surface nodes
bId = O.bodies.append(levelSetBody(clump=clump, spacing=spac * 1.e-6, nSurfNodes=nN, dynamic=True, center=(0, 0, zC), orientation=clumpB.state.ori))
lsClump = O.bodies[bId]  # lsClump body

# Wall z=0 ground:
bId = O.bodies.append(wall(0, 2))
floor = O.bodies[bId]  # floor body

###########
# Engines #
###########

knVal = 6.e4
ksVal = 0.5 * knVal
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_LevelSet_Aabb(), Bo1_Wall_Aabb()], verletDist=0),
        InteractionLoop(
                [Ig2_Wall_LevelSet_MultiScGeom()], [Ip2_FrictMat_FrictMat_MultiFrictPhys(kn=knVal, ks=ksVal)],
                [Law2_MultiScGeom_MultiFrictPhys_CundallStrack(sphericalBodies=False)]
        ),
        NewtonIntegrator(gravity=Vector3(0, 0, -9.8)),
        PyRunner(command='saveData()', iterPeriod=1, initRun=True),
        VTKRecorder(recorders=["lsBodies"], iterPeriod=100, fileName='multiScGeom', initRun=True)
]


def saveData():
	if O.interactions.has(floor.id, lsClump.id, True):
		cont = O.interactions[floor.id, lsClump.id]
		nC = len(cont.geom.contacts)
	else:
		nC = 0
	plot.addData(
	        it=O.iter,
	        nC=nC,
	        fFromForces=O.forces.f(
	                lsClump.id
	        )  # We will automatically get fFromForces_x, fFromForces_y, fFromForces_z in plot.data. Note that sync = True would be required with a true Clump
	        ,
	        pos=lsClump.state.pos
	)


O.dt = 2.97863676629746e-05  # 0.4*(mass/ kn)**0.5

#######
# Run #
#######

O.run(2001, True)  # 2001 should take us to equilibrium

########################
# Some post-processing #
########################

# graph
plot.plots = {'it': ('nC'), ' it': ('fFromForces_x', 'fFromForces_y', 'fFromForces_z'), 'it ': ('pos_x', 'pos_y'), ' it ': 'pos_z'}
plot.plot()

# back-verification of a correct force applied onto the LSclump
contactingPoints = O.interactions[floor.id, lsClump.id].geom.contacts  # all kinematic (ScGeom) items for the floor - lsClump interaction
forces = [
        (knVal * 2 * contPt.contactPoint[2] * Vector3.UnitZ) for contPt in contactingPoints
]  # sum kn un, while applying a *2 to contPt altitude because contactPoint is in the middle of overlap
fOnTop = numpy.sum(forces, axis=0)  # summing all elements encountered when walking along axis 0 (the only one possible to encounter Vector3)
print('\nWe should have equality of the following z-forces (modulo sign)', fOnTop[2], O.forces.f(5)[2])
print('Indeed, one + another =', fOnTop[2] + O.forces.f(5)[2], '\n')
