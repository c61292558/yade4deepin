// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2004 Janek Kozicki <cosurgi@berlios.de>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
// © 2006 Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>

#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((ScGeom)(ScGeom6D)(ChCylGeom6D));

Vector3r& ScGeom::rotate(Vector3r& shearForce) const
{
	// approximated rotations
	shearForce -= shearForce.cross(orthonormal_axis);
	shearForce -= shearForce.cross(twist_axis);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	//shearForce -= normal.dot(shearForce)*normal;
	return shearForce;
}


Vector3r& ScGeom::rotateNonSpherical(Vector3r& shearForce) const
{ //FIXME
	// approximated rotations
	shearForce -= shearForce.cross(orthonormal_axis);
	//shearForce -= shearForce.cross(twist_axis);
	//NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	//shearForce -= normal.dot(shearForce)*normal;
	if (math::isnan(shearForce.norm())) { std::cout << "orthonormal_axis: " << orthonormal_axis << ", normal: " << normal << endl; }
	return shearForce;
}

//!Precompute data needed for rotating tangent vectors attached to the interaction
void ScGeom::precompute(
        const State&                   rbp1,
        const State&                   rbp2,
        const Scene*                   scene,
        const shared_ptr<Interaction>& c,
        const Vector3r&                currentNormal,
        bool                           isNew,
        const Vector3r&                shift2,
        bool                           avoidGranularRatcheting)
{
	if (!isNew) {
		orthonormal_axis = normal.cross(currentNormal);
		Real angle       = scene->dt * 0.5 * normal.dot(rbp1.angVel + rbp2.angVel);
		twist_axis       = angle * normal;
	} else
		twist_axis = orthonormal_axis = Vector3r::Zero();
	//Update contact normal
	normal = currentNormal;
	//Precompute shear increment
	Vector3r relativeVelocity = getIncidentVel(
	        &rbp1, &rbp2, scene->dt, shift2, scene->isPeriodic ? scene->cell->intrShiftVel(c->cellDist) : Vector3r::Zero(), avoidGranularRatcheting);
	//keep the shear part only
	relativeVelocity = relativeVelocity - normal.dot(relativeVelocity) * normal;
	shearInc         = relativeVelocity * scene->dt;
}

void ScGeom::doIg2Work(
        Vector3r                       ctctPt,
        Real                           un,
        Real                           rad1,
        Real                           rad2,
        const State&                   rbp1,
        const State&                   rbp2,
        const Scene*                   scene,
        const shared_ptr<Interaction>& c,
        const Vector3r&                currentNormal,
        const Vector3r&                shift2,
        bool                           newScGeom,
        bool                           avoidGranularRatcheting)
{
	// goal is to avoid duplicating the definition of ScGeom attributes and execution of its ::precompute.
	contactPoint     = ctctPt;
	penetrationDepth = un;
	// NB radius1, radius2: those are useful for
	// 1. contact kinematics description if and only if avoidGranularRatcheting
	// 2. applying contact forces in C-S Law2, if sphericalBodies
	// 3. GSTS time step determination with respect to rotational stiffnesses
	// 4. and also, after being coined as refR1, refR2, for contact stiffness expression in FrictPhys/FrictMat
	radius1 = rad1;
	radius2 = rad2;
	precompute(rbp1, rbp2, scene, c, currentNormal, newScGeom, shift2, avoidGranularRatcheting);
	// precompute will take care of
	// 1. preparing the rotation of shearForce to the new tangent plane (done later, in Law2) defining these orthonormal_axis and twist_axis, if newScGeom = false
	// 2. updating geomPtr->normal (previous value) to currentNormal
	// 3. computing the relative velocity at contact, through getIncidentVel(avoidGranularRatcheting), using now-defined contactPoint
	//	Comparing with Ig2_Sphere_Sphere_ScGeom.cpp, I think everything is here..
}

Vector3r
ScGeom::getIncidentVel(const State* rbp1, const State* rbp2, Real /*dt*/, const Vector3r& shift2, const Vector3r& shiftVel, bool avoidGranularRatcheting) const
{
	if (avoidGranularRatcheting) {
		/* B.C. Comment :
		Short explanation of what we want to avoid :
		Numerical ratcheting is best understood considering a small elastic cycle at a contact between two grains : assuming b1 is fixed, impose this displacement to b2 :
		1. translation "dx" in the normal direction
		2. rotation "a"
		3. translation "-dx" (back to initial position)
		4. rotation "-a" (back to initial orientation)
		If the branch vector used to define the relative shear in rotation×branch is not constant (typically if it is defined from the vector center→contactPoint), then the shear displacement at the end of this cycle is not zero: rotations *a* and *-a* are multiplied by branches of different lengths.
		It results in a finite contact force at the end of the cycle even though the positions and orientations are unchanged, in total contradiction with the elastic nature of the problem. It could also be seen as an *inconsistent energy creation or loss*. Given that DEM simulations tend to generate oscillations around equilibrium (damped mass-spring), it can have a significant impact on the evolution of the packings, resulting for instance in slow creep in iterations under constant load.
		The solution adopted here to avoid ratcheting is as proposed by McNamara and co-workers.
		They analyzed the ratcheting problem in detail - even though they comment on the basis of a cycle that differs from the one shown above. One will find interesting discussions in e.g. DOI 10.1103/PhysRevE.77.031304, even though solution it suggests is not fully applied here (equations of motion are not incorporating alpha, in contradiction with what is suggested by McNamara et al.).
		 */
		// In case one of two bodies is not actually Sphere-shaped, this will give an erroneous value of relative velocity...
		Real     alpha            = (radius1 + radius2) / (radius1 + radius2 - penetrationDepth);
		Vector3r relativeVelocity = (rbp2->vel - rbp1->vel) * alpha + rbp2->angVel.cross(-radius2 * normal) - rbp1->angVel.cross(radius1 * normal);
		relativeVelocity += alpha * shiftVel;
		return relativeVelocity;
	} else {
		// Here, the following is always correct (even for non Sphere-shaped bodies)
		Vector3r c1x              = (contactPoint - rbp1->pos);
		Vector3r c2x              = (contactPoint - (rbp2->pos + shift2));
		Vector3r relativeVelocity = (rbp2->vel + rbp2->angVel.cross(c2x)) - (rbp1->vel + rbp1->angVel.cross(c1x));
		relativeVelocity += shiftVel;
		return relativeVelocity;
	}
}

Vector3r ScGeom::getIncidentVel(const State* rbp1, const State* rbp2, Real dt, bool avoidGranularRatcheting) const
{
	//Just pass null shift to the periodic version
	return getIncidentVel(rbp1, rbp2, dt, Vector3r::Zero(), Vector3r::Zero(), avoidGranularRatcheting);
}

Vector3r ScGeom::getIncidentVel_py(shared_ptr<Interaction> i, bool avoidGranularRatcheting) const
{
	if (i->geom.get() != this) throw invalid_argument("ScGeom object is not the same as Interaction.geom.");
	Scene* scene = Omega::instance().getScene().get();
	return getIncidentVel(
	        Body::byId(i->getId1(), scene)->state.get(),
	        Body::byId(i->getId2(), scene)->state.get(),
	        scene->dt,
	        scene->isPeriodic ? scene->cell->intrShiftPos(i->cellDist) : Vector3r::Zero(), // shift2
	        scene->isPeriodic ? scene->cell->intrShiftVel(i->cellDist) : Vector3r::Zero(), // shiftVel
	        avoidGranularRatcheting);
}

Vector3r ScGeom::getRelAngVel(const State* rbp1, const State* rbp2, Real /*dt*/) const
{
	Vector3r relAngVel = (rbp2->angVel - rbp1->angVel);
	return relAngVel;
}

Vector3r ScGeom::getRelAngVel_py(shared_ptr<Interaction> i) const
{
	if (i->geom.get() != this) throw invalid_argument("ScGeom object is not the same as Interaction.geom.");
	Scene* scene = Omega::instance().getScene().get();
	return getRelAngVel(Body::byId(i->getId1(), scene)->state.get(), Body::byId(i->getId2(), scene)->state.get(), scene->dt);
}

//!Precompute relative rotations (and precompute ScGeom3D)
void ScGeom6D::precomputeRotations(const State& rbp1, const State& rbp2, bool isNew, bool creep)
{
	if (isNew) {
		initRotations(rbp1, rbp2);
	} else {
		Quaternionr delta((rbp1.ori * (initialOrientation1.conjugate())) * (initialOrientation2 * (rbp2.ori.conjugate())));
		delta.normalize();
		if (creep) delta = delta * twistCreep;
		AngleAxisr aa(delta); // axis of rotation - this is the Moment direction UNIT vector; // angle represents the power of resistant ELASTIC moment
		                      //Eigen::AngleAxisr(q) returns nan's when q close to identity, next tline fixes the pb.
// add -DYADE_SCGEOM_DEBUG to CXXFLAGS to enable this piece or just do
// #define YADE_SCGEOM_DEBUG //(but do not commit with that enabled in the code)
#ifdef YADE_SCGEOM_DEBUG
		if (math::isnan(aa.angle())) {
			cerr << "NaN angle found in angleAxisr(q), for quaternion " << delta << ", after quaternion product" << endl;
			cerr << "rbp1.ori * (initialOrientation1.conjugate())) * (initialOrientation2 * (rbp2.ori.conjugate()) with quaternions :" << endl;
			cerr << rbp1.ori << " * " << initialOrientation1 << " * " << initialOrientation2 << " * " << rbp2.ori << endl
			     << " and sub-products :" << endl
			     << rbp1.ori * (initialOrientation1.conjugate()) << " * " << initialOrientation2 * (rbp2.ori.conjugate()) << endl;
			cerr << "q.w (before normalization) " << delta.w();
			delta.normalize();
			cerr << "q.w (after) " << delta.w() << endl;
			AngleAxisr bb(delta);
			cerr << delta << " " << bb.angle() << endl;
		}
#else
		if (math::isnan(aa.angle())) aa.angle() = 0;
#endif
		if (aa.angle() > Mathr::PI) aa.angle() -= Mathr::TWO_PI; // angle is between 0 and 2*pi, but should be between -pi and pi
		twist   = (aa.angle() * aa.axis().dot(normal));
		bending = Vector3r(aa.angle() * aa.axis() - twist * normal);
	}
}

void ScGeom6D::initRotations(const State& state1, const State& state2)
{
	initialOrientation1 = state1.ori;
	initialOrientation2 = state2.ori;
	twist               = 0;
	bending             = Vector3r::Zero();
	twistCreep          = Quaternionr(1.0, 0.0, 0.0, 0.0);
}

} // namespace yade
