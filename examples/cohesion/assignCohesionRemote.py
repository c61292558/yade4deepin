"""
This exemple originated from a discussion in yade-dev/answers on Gitlab [1].
It shows how to create cohesive bonds between distant particles when the initial configuration
is assumed force-free.
[1] https://gitlab.com/yade-dev/answers/-/issues/4022
"""



from yade import pack, plot
import math

#Define materials
cohMat = CohFrictMat(
    young=1e7,
    poisson=0.3,
    frictionAngle=radians(30),
    density=2600,
    isCohesive=True,
    normalCohesion=1e5,
    shearCohesion=1e5
)

#Create spheres
rad = 0.05
#desiredEquilibriumDistance = 0.02 #Desired equilibrium distance between sphere centers
initialSeparation = 0.01  # Initial separation beyond equilibrium
s1 = sphere([0, 0, 0], radius=rad, material=cohMat, fixed=True)
s2 = sphere([2*rad+initialSeparation, 0, 0], radius=rad, material=cohMat)
s3 = sphere([3.5*rad+2*initialSeparation, 0, 0], radius=0.5*rad, material=cohMat)

#Add initial velocity to s2
s2.state.vel = Vector3(0.1, 0, 0)  # 0.1 m/s in positive x direction

#Add spheres to the scene
O.bodies.append((s1, s2, s3))

#Define engines
O.engines = [
    ForceResetter(),
#we enlarge detection distance for the collider to detect all neighbors within a certain range
#it is later reverted for performance                                                    reasons.
    InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=1.5)],label="collider"),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom6D(interactionDetectionFactor=1.5,label="ig2")],
        [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
            setCohesionNow=False,
            setCohesionOnNewContacts=False,
            label = "physFunctor"
        )],
        [Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
    ),
    NewtonIntegrator(damping=0),
    PyRunner(command='addPlotData()', iterPeriod=1)
]

#Function to add plot data
def addPlotData():
    plot.addData(n=O.iter, fn=O.forces.f(0)[0], un=s2.state.pos[0]-s1.state.pos[0]-2*rad)

#Set up plot
plot.plots = {'n': ('fn', None,'un')}
plot.plot()

#Set time step
O.dt = 5e-5

#detect all interactions within 1.5 radius distance
collider.__call__()

#loop over the interactions found by collider
#they are all virtual, so we need "all()" to see them iterable
for i in O.interactions.all():
#you may want additional conditions here depending on the use case
    if (True):
        i=createInteraction(i.id1,i.id2) # now it is real, not virtual
        physFunctor.setCohesion(i,cohesive=True,resetDisp=True)


## Run simulation

#No need to detect distant neighbors any more, the collider will be faster like this 
collider.aabbEnlargeFactor = 1
ig2.interactionDetectionFactor = 1

O.run(2000, True)
