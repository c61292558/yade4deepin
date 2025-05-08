These scripts illustrate the use of the [LevelSet](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.LevelSet) class for shape description.

Please note in particular:

- `discharge.py` for a LevelSet-based simulation of a discharge of 1000 superquadric particles, as discussed in Section 6 of [[Duriez2021b](https://www.sciencedirect.com/science/article/pii/S0098300421002247)]
- `levelSetBody.py` for a particle-scale (1 or 2 bodies) simulation in order to illustrate the various syntaxes available for defining a LevelSet-shaped body in YADE and the expected `O.engines`

It is furthermore provided 3 utilities functions in separate files:

- `pvVisu.py` for a Python-automated vizualization of LevelSet-shaped YADE bodies in Paraview (after using `VTKRecorder` in the simulation, see `levelSetBody.py` for an example)
- `ramFP.py` for measuring current RAM usage on your machine, as this could be a limiting factor when using LevelSet
- `spaceRot.py` for defining "random" orientations in 3D space