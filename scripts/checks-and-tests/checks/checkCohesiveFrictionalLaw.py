from numpy import array
#from yade import plot   ## uncomment plotting lines for debugging

O.periodic = True
O.cell.hSize = Matrix3(10, 0, 0, 0, 10, 0, 0, 0, 10)

results = []


def getState():
	state = [float(O.interactions.has(0, 1) and O.interactions[0, 1].isReal)]
	if state[0]:
		i = O.interactions[0, 1]
		for x in [
		        int(i.phys.cohesionBroken),
		        i.phys.normalForce.norm(),
		        i.phys.shearForce.norm(),
		        i.phys.moment_bending.norm(),
		        i.phys.moment_twist.norm()
		]:
			state.append(x)
	else:
		state = state + [
		        0, 0, 0, 0, 0
		]  # not all numpy versions will handle nicely 2D arrays with heterogeneous lists, make sure all states have 6 floats
	results.append(numpy.array(state))
	print(state)


expectedResults = [
        array([1.0, 0, 0.0, 0.0, 0.0, 0.0]),
        array([1.0, 0, 399.5999999999995, 0.0, 0.0, 0.0]),
        array([1.0, 0, 899.0999999999989, 0.0, 0.0, 0.0]),
        array([0.0, 0, 0, 0, 0, 0]),
        array([0.0, 0, 0, 0, 0, 0]),
        array([1.0, 1, 98.89999999999866, 0.0, 0.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 0.0, 10.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 0.0, 60.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 0.0, 0.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 0.0, 0.0, 98.89999999999866]),
        array([1.0, 1, 98.89999999999866, 0.0, 50.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 0.0, 98.89999999999866, 0.0]),
        array([1.0, 1, 98.89999999999866, 30.0, 0.0, 0.0]),
        array([1.0, 1, 98.89999999999866, 54.02931624555014, 0.0, 0.0])
]

O.engines = [  # define engines, main functions for simulation
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()],
      label='ISCollider'),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom6D()],
                [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=True,label="ip2")],
                [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm=True,always_use_moment_law=True)]
        ),
        #PyRunner(command='plot.addData(ii=O.iter+ref,x=O.bodies[1].state.pos[0]-7,Fx=O.forces.f(0)[0],Mz=O.forces.t(0)[2])', iterPeriod=1,initRun=True),
        NewtonIntegrator(label='newton')
]

O.dt = 1
ref = 0
O.materials.append(CohFrictMat(young=1e3, normalCohesion=1e3, shearCohesion=0.5e3, etaRoll=1, etaTwist=1, momentRotationLaw=True, fragile=True))
O.bodies.append(sphere((4, 5, 5), 1, fixed=True))
O.bodies.append(sphere((7, 5, 5), 1, fixed=True))

i = createInteraction(0, 1)
ip2.setCohesionOnNewContacts = False
i.phys.unp = -(O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos).norm() + O.bodies[i.id2].shape.radius + O.bodies[i.id1].shape.radius

O.saveTmp()
#plot.plots={'ii':('Fx',None,'Mz')}
#plot.plot()


def init():
	global ref
	ref += O.iter
	O.reload()
	O.step()  # before any motion so one data point gives the origin when plotting


init()
print("Pull contact until it breaks")
getState()
O.bodies[1].state.vel = (0.0999, 0, 0)
O.run(5, True)
getState()
O.run(5, True)
getState()
O.run(2, True)
getState()
print("Inverse velocity to create new non-cohesive contact")
O.bodies[1].state.vel = -O.bodies[1].state.vel
O.run(23, True)
getState()
O.bodies[1].state.vel = (0, 0, 0)
O.run(1, True)
getState()
O.saveTmp()
print("Apply pure bending")
O.bodies[1].state.angVel = (0, 0, 0.01)
O.bodies[0].state.angVel = (0, 0, -0.01)
O.run(1, True)
getState()
O.run(5, True)
getState()
print("Un-bend")
O.bodies[1].state.angVel = (0, 0, -0.01)
O.bodies[0].state.angVel = (0, 0, 0.01)
O.run(6, True)
getState()
O.bodies[0].state.angVel = O.bodies[1].state.angVel = (0, 0, 0)
O.step()
print("Twist")
O.bodies[0].state.angVel = (0.02, 0, 0)
O.run(20, True)
getState()
print("Break remote interaction by bending")
init()
O.bodies[1].state.angVel = (0, 0, 0.01)
O.bodies[0].state.angVel = (0, 0, -0.01)
O.run(5, True)
getState()
O.run(10, True)
getState()
print("Shear compressed contact with residual friction")
init()
O.bodies[1].state.angVel = (0, 0, 0.02)
O.bodies[0].state.angVel = (0, 0, 0.02)
O.run(3, True)
getState()
O.run(10, True)
getState()

difference = numpy.linalg.norm(numpy.array(results) - numpy.array(expectedResults))
if difference > 1e-8:
	raise YadeCheckError('Results of CohesiveFrictional law have changed (difference =' + str(numpy.array(results) - numpy.array(expectedResults)) + ')')
else:
	print('CohesiveFrictional model passed with difference =', difference)
