# -*- encoding=utf-8 -*-
# create a few clump configurations by hand
from yade import pack

c1 = pack.SpherePack([((0, 0, 0), .5), ((.5, 0, 0), .5), ((0, .5, 0), .3)])
c2 = pack.SpherePack([((0, 0, 0), .5), ((.7, 0, 0), .3), ((.9, 0, 0), .2)])
sp = pack.SpherePack()
print('Generated # of clumps:', sp.makeClumpCloud((0, 0, 0), (15, 15, 15), [c1, c2], periodic=False))
sp.toSimulation()

O.bodies.append(wall(position=0, axis=2))

O.engines = [
        #SubdomainBalancer(),
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        NewtonIntegrator(damping=.4, gravity=(0, 0, -100))
]
O.dt = .7 * PWaveTimeStep()
O.saveTmp()
O.step()
