import os
# from yadeimport import *
from yade import mpy as mp

### example to test the hydrodynamic forces and correct OF-YADE coupling  ###

compFricCoef = 0.0 # initial contact friction during the confining phase
finalFricCoef = 0.47 # contact friction during the deviatoric loading
wallFricCoef = 0.19 # contact friction during the deviatoric loading
densitySpheres = 1500
g = 9.81

parallelYade=False #mpirun --allow-run-as-root -n 2 python3 scriptMPI.py , if False  python3 scriptMPI.py
numProcOF=2

young=6.3e7
poissonR=0.5
L_tank = 0.04
L = 0.01
H = 0.02
Radius=0.0005
porosity=0.6
restitCoef=0.91

density = 1000
NSTEPS = 10000000
saveVTK=1000
O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=densitySpheres, frictionAngle=compFricCoef, label='spheremat'))
O.materials.append(ViscElMat(en=restitCoef, et=1., young=young, poisson=poissonR, density=0, frictionAngle=wallFricCoef, label='walls'))

O.bodies.append(sphere((L/2., H*3/4., L/2.),Radius, material='spheremat'))

print ("Numer of particles:", len(O.bodies))


fluidCoupling = FoamCoupling()
fluidCoupling.couplingModeParallel = parallelYade
fluidCoupling.isGaussianInterp = True
#use pimpleFoamYade for gaussianInterp (only in serial mode)
sphereIDs = [b.id for b in O.bodies if type(b.shape) == Sphere]


'''The yade specific (icoFoamYade, pimpleFoamYade) OpenFOAM solver can be found in $FOAM_USER_APPBIN, (
# full path here, the scond argument, 2 is the number of FoamProcs. '''
# fluidCoupling.SetOpenFoamSolver(os.environ.get('FOAM_USER_APPBIN')+'/icoFoamYade', 2)
# it also work without path after sourcing OFoam's bashrc
fluidCoupling.SetOpenFoamSolver("pimpleFoamYade", numProcOF)


with open('data.txt', 'a') as f:
            f.write("Time" + " " +"velocity"  + " " +"position" +  "\n")
            f.close


def saveData():
    for b in O.bodies:
        if type(b.shape) == Sphere:
                with open('data.txt', 'a') as f:
                            f.write(str(O.time) + " "+ str(b.state.vel[1]) + " " +str(b.state.pos[1]) +  "\n")
                            f.close
    return



newton=NewtonIntegrator(gravity=(-0,-g,0),damping=0.4)

O.engines=[
 ForceResetter(),
 InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], label='collider', allowBiggerThanPeriod=True),
 #InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Wall_Aabb()]),
 InteractionLoop(
  [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],
  #[Ig2_Sphere_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor),Ig2_Box_Sphere_ScGeom6D(interactionDetectionFactor=enlFactor), Ig2_Wall_Sphere_ScGeom()],
  [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
  [Law2_ScGeom_ViscElPhys_Basic()],label="InteractionLoop"
 ),
        GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.7, timeStepUpdateInterval=100, parallelMode=parallelYade, label="ts"),
        fluidCoupling,  #to be called after timestepper
        newton,
        VTKRecorder(fileName='spheres/3d-vtk-', recorders=['spheres','boxes'], parallelMode=parallelYade, iterPeriod=saveVTK),
        PyRunner(iterPeriod=10,command='saveData()',dead=0,label="data")]


collider.verletDist = 0.0075
mp.YADE_TIMING = False
mp.FLUID_COUPLING = True
mp.VERBOSE_OUTPUT = False
mp.USE_CPP_INTERS = True
mp.ERASE_REMOTE_MASTER = True
mp.REALLOC_FREQUENCY = 0
mp.fluidBodies = sphereIDs
#mp.commSplit = True
mp.DOMAIN_DECOMPOSITION = True


O.dt = 1e-5
O.dynDt = False

# dataFProfile.dead=0

mp.mpirun(NSTEPS)
mp.mprint("RUN FINISH")
fluidCoupling.killMPI()

exit()
