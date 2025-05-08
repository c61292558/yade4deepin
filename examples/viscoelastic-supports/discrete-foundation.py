import numpy as np
from yade import plot, qt

# PARAMS
box_x = 0.05
box_y = 0.05
box_z = 0.05
box_A = box_x * box_y
box_grid_size = 40

k1 = 10e5 * box_A
k2 = 10e5 * box_A
c1 = 100e4 * box_A
c2 = 100e2 * box_A

g = 10  # gravity constant
# spheres
r = 0.75
H = 1.0
# BODIES
# Boxes
box_ids = []
for i in range(box_grid_size):
	for j in range(box_grid_size):
		b_id = O.bodies.append(box(center=(i * box_x, j * box_y, 0), extents=(box_x / 2, box_y / 2, box_z / 2), dynamic=True, color=(0, 0, 1)))
		box_ids += [b_id]
		O.bodies[b_id].state.blockedDOFs = 'xyXYZ'

		# add 'geostatic stress' to remove effect of gravity on the foundation (otherwise we would have to let the foundation settle first0
		mass = O.bodies[b_id].state.mass
		O.forces.setPermF(b_id, (0, 0, mass * g))

# sphere
sph_id = O.bodies.append(sphere(center=((1.01 + box_grid_size * box_x) / 2, box_grid_size * box_y / 2, box_z / 2 + H), radius=r))

# ENGINES

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()],
        ),
        VESupportEngine(bIds=box_ids, k1=k1, k2=k2, c1=c1, c2=c2),
        NewtonIntegrator(damping=0.0, gravity=(0, 0, -g), label='newton'),
        PyRunner(command="colorize_by_disp(ids= box_ids)", iterPeriod=25),
        PyRunner(command="plot.addData(t =O.time, disp =O.bodies[sph_id].state.pos[2])", iterPeriod=100)
]

# SIMULATION CONTROL

O.dt = 1.e-4
#O.dt = 0.8*PWaveTimeStep()


def colorize_by_disp(ids, min_scalar=0, max_scalar=0.15):
	for b_id in ids:
		displ = O.bodies[b_id].state.displ().norm()
		color = scalarOnColorScale(displ, min_scalar, max_scalar)
		O.bodies[b_id].shape.color = color


plot.plots = {'t': 'disp'}

plot.plot()

#view
v = qt.View()
v.upVector = (0, -0.25, 1)
v.showEntireScene()

O.run(50000)
