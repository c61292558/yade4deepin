from yade import plot, qt

# PARAMS
box_x = 0.1
box_y = 0.1
box_z = 0.0001

box_A = box_x * box_y

k1 = 10e5 * box_A
k2 = 10e5 * box_A
c1 = -1 * box_A
c2 = 100e3 * box_A

g = 10
# spheres
r = 0.1
H = 1  # packing height

# BODIES
# Box
b_id = O.bodies.append(box(center=(0, 0, -box_z / 2), extents=(box_x / 2, box_y / 2, box_z / 2), dynamic=True, color=(0, 0, 1)))
O.bodies[b_id].state.blockedDOFs = 'xyXYZ'

# Sphere
s = O.bodies.append(sphere(center=(0, 0, H), radius=r))

# Material
O.materials[0].young = 1e5

# ENGINES

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        VESupportEngine(bIds=[b_id], k1=k1, k2=k2, c1=c1, c2=c2),
        NewtonIntegrator(damping=0.0, gravity=(0, 0, -g), label='newton'),
        PyRunner(command='addPlotData(stop_time = 7)', firstIterRun=1, iterPeriod=100),
]

# SIMULATION CONTROL

O.dt = 1.e-6


def addPlotData(stop_time):
	f1 = O.forces.f(s)[2]
	h1 = O.bodies[s].state.pos[2] - r
	h2 = O.bodies[b_id].state.pos[2] + box_z / 2

	plot.addData(t=O.time, f1=f1, h1=h1, h2=h2)

	if O.time > stop_time:
		O.pause()
		plot.saveDataTxt('data_example_03.txt')


plot.plots = {'t': ('f1', None, 'h1', 'h2')}
plot.plot()
#view
v = qt.View()
v.upVector = (0, 0, 1)
v.showEntireScene()

#O.run()
