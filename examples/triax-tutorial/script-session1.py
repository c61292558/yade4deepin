# -*- encoding=utf-8 -*-
#*************************************************************************
#  Copyright (C) 2010 by Bruno Chareyre                                  *
#  bruno.chareyre_at_grenoble-inp.fr                                     *
#                                                                        *
#  This program is free software; it is licensed under the terms of the  *
#  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/

## This script details the simulation of a triaxial test on sphere packings using Yade
## See the associated pdf file for detailed exercises
## the algorithms presented here have been used in published papers, namely:
## * Chareyre et al. 2002 (http://www.geosyntheticssociety.org/Resources/Archive/GI/src/V9I2/GI-V9-N2-Paper1.pdf)
## * Chareyre and Villard 2005 (https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
## * Scholtès et al. 2009 (http://dx.doi.org/10.1016/j.ijengsci.2008.07.002)
## * Tong et al.2012 (http://dx.doi.org/10.2516/ogst/2012032)
##
## Most of the ideas were actually developped during my PhD.
## If you want to know more on micro-macro relations evaluated by triaxial simulations
## AND if you can read some french, it is here: http://tel.archives-ouvertes.fr/docs/00/48/68/07/PDF/Thesis.pdf

from yade import pack

############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################

# The following 5 lines will be used later for batch execution
nRead = readParamsFromTable(
        num_spheres=1000,  # number of spheres
        compFricDegree=30,  # contact friction during the confining phase
        key='_triax_base_',  # put you simulation's name here
        unknownOk=True
)
from yade.params import table

num_spheres = table.num_spheres  # number of spheres
key = table.key
targetPorosity = 0.43  #the porosity we want for the packing
compFricDegree = table.compFricDegree  # initial contact friction during the confining phase (will be decreased during the REFD compaction process)
finalFricDegree = 30  # contact friction during the deviatoric loading
rate = -0.02  # loading rate (strain rate)
damp = 0.2  # damping coefficient
stabilityThreshold = 0.01  # we test unbalancedForce against this value in different loops (see below)
young = 5e6  # contact stiffness
mn, mx = Vector3(0, 0, 0), Vector3(1, 1, 1)  # corners of the initial packing

## create materials for spheres and plates
O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=radians(compFricDegree), density=2600, label='spheres'))
O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=0, density=0, label='walls'))

## create walls around the packing
walls = aabbWalls([mn, mx], thickness=0, material='walls')
wallIds = O.bodies.append(walls)

## use a SpherePack object to generate a random loose particles packing
sp = pack.SpherePack()

clumps = False  #turn this true for the same example with clumps
if clumps:
	## approximate mean rad of the futur dense packing for latter use
	volume = (mx[0] - mn[0]) * (mx[1] - mn[1]) * (mx[2] - mn[2])
	mean_rad = pow(0.09 * volume / num_spheres, 0.3333)
	## define a unique clump type (we could have many, see clumpCloud documentation)
	c1 = pack.SpherePack([((-0.2 * mean_rad, 0, 0), 0.5 * mean_rad), ((0.2 * mean_rad, 0, 0), 0.5 * mean_rad)])
	## generate positions and input them in the simulation
	sp.makeClumpCloud(mn, mx, [c1], periodic=False)
	sp.toSimulation(material='spheres')
	O.bodies.updateClumpProperties()  #get more accurate clump masses/volumes/inertia
else:
	sp.makeCloud(mn, mx, -1, 0.3333, num_spheres, False, 0.95, seed=1)  #"seed" make the "random" generation always the same
	O.bodies.append([sphere(center, rad, material='spheres') for center, rad in sp])
	#or alternatively (higher level function doing exactly the same):
	#sp.toSimulation(material='spheres')

############################
###   DEFINING ENGINES   ###
############################

triax = TriaxialStressController(
        ## TriaxialStressController will be used to control stress and strain. It controls particles size and plates positions.
        ## this control of boundary conditions was used for instance in http://dx.doi.org/10.1016/j.ijengsci.2008.07.002
        maxMultiplier=1. + 2e4 / young,  # spheres growing factor (fast growth)
        finalMaxMultiplier=1. + 2e3 / young,  # spheres growing factor (slow growth)
        thickness=0,
        ## switch stress/strain control using a bitmask. What is a bitmask, huh?!
        ## Say x=1 if stess is controlled on x, else x=0. Same for for y and z, which are 1 or 0.
        ## Then an integer uniquely defining the combination of all these tests is: mask = x*1 + y*2 + z*4
        ## to put it differently, the mask is the integer whose binary representation is xyz, i.e.
        ## "100" (1) means "x", "110" (3) means "x and y", "111" (7) means "x and y and z", etc.
        stressMask=7,
        internalCompaction=True,  # If true the confining pressure is generated by growing particles
)

newton = NewtonIntegrator(damping=damp)

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        ## We will use the global stiffness of each body to determine an optimal timestep (see https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
        GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=100, timestepSafetyCoefficient=0.8),
        triax,
        TriaxialStateRecorder(iterPeriod=100, file='WallStresses' + table.key),
        newton
]

if nRead == 0:
	yade.qt.Controller(), yade.qt.View()

## UNCOMMENT THE FOLLOWING SECTIONS ONE BY ONE
## DEPENDING ON YOUR EDITOR, IT COULD BE DONE
## BY SELECTING THE CODE BLOCKS BETWEEN THE SUBTITLES
## AND PRESSING CTRL+SHIFT+D

#######################################
###   APPLYING CONFINING PRESSURE   ###
#######################################

#the value of (isotropic) confining stress defines the target stress to be applied in all three directions
triax.goal1 = triax.goal2 = triax.goal3 = -10000

#while 1:
#O.run(1000, True)
##the global unbalanced force on dynamic bodies, thus excluding boundaries, which are not at equilibrium
#unb=unbalancedForce()
#print 'unbalanced force:',unb,' mean stress: ',triax.meanStress
#if unb<stabilityThreshold and abs(-10000-triax.meanStress)/10000<0.001:
#break

#O.save('confinedState'+key+'.yade.gz')
#print "###      Isotropic state saved      ###"

###################################################
###   REACHING A SPECIFIED POROSITY PRECISELY   ###
###################################################

### We will reach a prescribed value of porosity with the REFD algorithm
### (see http://dx.doi.org/10.2516/ogst/2012032 and
### http://www.geosyntheticssociety.org/Resources/Archive/GI/src/V9I2/GI-V9-N2-Paper1.pdf)

#import sys #this is only for the flush() below
#while triax.porosity>targetPorosity:
## we decrease friction value and apply it to all the bodies and contacts
#compFricDegree = 0.95*compFricDegree
#setContactFriction(radians(compFricDegree))
#print "\r Friction: ",compFricDegree," porosity:",triax.porosity,
#sys.stdout.flush()
## while we run steps, triax will tend to grow particles as the packing
## keeps shrinking as a consequence of decreasing friction. Consequently
## porosity will decrease
#O.run(500,1)

#O.save('compactedState'+key+'.yade.gz')
#print "###    Compacted state saved      ###"

##############################
###   DEVIATORIC LOADING   ###
##############################

##We move to deviatoric loading, let us turn internal compaction off to keep particles sizes constant
#triax.internalCompaction=False

## Change contact friction (remember that decreasing it would generate instantaneous instabilities)
#setContactFriction(radians(finalFricDegree))

##set stress control on x and z, we will impose strain rate on y
#triax.stressMask = 5
##now goal2 is the target strain rate
#triax.goal2=rate
## we define the lateral stresses during the test, here the same 10kPa as for the initial confinement.
#triax.goal1=-10000
#triax.goal3=-10000

##we can change damping here. What is the effect in your opinion?
#newton.damping=0.1

##Save temporary state in live memory. This state will be reloaded from the interface with the "reload" button.
#O.saveTmp()

#####################################################
###    Example of how to record and plot data     ###
#####################################################

#from yade import plot

### a function saving variables
#def history():
#plot.addData(e11=-triax.strain[0], e22=-triax.strain[1], e33=-triax.strain[2],
#ev=-triax.strain[0]-triax.strain[1]-triax.strain[2],
#s11=-triax.stress(triax.wall_right_id)[0],
#s22=-triax.stress(triax.wall_top_id)[1],
#s33=-triax.stress(triax.wall_front_id)[2],
#i=O.iter)

#if 1:
## include a periodic engine calling that function in the simulation loop
#O.engines=O.engines[0:5]+[PyRunner(iterPeriod=20,command='history()',label='recorder')]+O.engines[5:7]
##O.engines.insert(4,PyRunner(iterPeriod=20,command='history()',label='recorder'))
#else:
## With the line above, we are recording some variables twice. We could in fact replace the previous
## TriaxialRecorder
## by our periodic engine. Uncomment the following line:
#O.engines[4]=PyRunner(iterPeriod=20,command='history()',label='recorder')

#O.run(100,True)

### declare what is to plot. "None" is for separating y and y2 axis
#plot.plots={'i':('e11','e22','e33',None,'s11','s22','s33')}
### the traditional triaxial curves would be more like this:
##plot.plots={'e22':('s11','s22','s33',None,'ev')}

## display on the screen (doesn't work on VMware image it seems)
#plot.plot()

#####  PLAY THE SIMULATION HERE WITH "PLAY" BUTTON OR WITH THE COMMAND O.run(N)  #####

## In that case we can still save the data to a text file at the the end of the simulation, with:
#plot.saveDataTxt('results'+key)
##or even generate a script for gnuplot. Open another terminal and type  "gnuplot plotScriptKEY.gnuplot:
#plot.saveGnuplot('plotScript'+key)
