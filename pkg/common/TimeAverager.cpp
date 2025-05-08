// 2024 © Benjamin Dedieu <benjamin.dedieu@proton.me>

#include "TimeAverager.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Scene.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/DemXDofGeom.hpp>
#include <boost/unordered_map.hpp>

// The TimeAverager class is used to average data over time for specific sphere identified by ids.
// See description in TimeAverager.hpp

namespace yade {

YADE_PLUGIN((TimeAverager));

// PUBLIC METHODS

Vector3r         TimeAverager::getPos(Body::id_t id) const { return getValueFromMap(pos, id); }
Vector3r         TimeAverager::getVel(Body::id_t id) const { return getValueFromMap(vel, id); }
Vector3r         TimeAverager::getAngVel(Body::id_t id) const { return getValueFromMap(angVel, id); }
Vector3r         TimeAverager::getForce(Body::id_t id) const { return getValueFromMap(force, id); }
Vector3r         TimeAverager::getTorque(Body::id_t id) const { return getValueFromMap(torque, id); }
Real             TimeAverager::getNbContact(Body::id_t id) const { return getValueFromMap(nbContact, id); }
Vector3r         TimeAverager::getContactForce(Body::id_t id) const { return getValueFromMap(contactForce, id); }
Vector3r         TimeAverager::getContactTorque(Body::id_t id) const { return getValueFromMap(contactTorque, id); }
vector<Vector3r> TimeAverager::getContactForceField(Body::id_t id) const
{
	if (computeContactForceField) { return getValueFromMap(contactForceField, id); }
	throw std::runtime_error("No value to retrieve for contactForceField since computeContactForceField is false.");
}

void TimeAverager::initialization()
{
	scene->forces.sync();
	FOREACH(Body::id_t id, ids)
	{
		const shared_ptr<Body>& b = Body::byId(id, scene);
		if (!b) continue;
		pos[id]           = b->state->pos;
		vel[id]           = b->state->vel;
		angVel[id]        = b->state->angVel;
		force[id]         = scene->forces.getForce(id);
		torque[id]        = scene->forces.getTorque(id);
		nbContact[id]     = getInstantNbContact(b);
		contactForce[id]  = getInstantContactForce(b);
		contactTorque[id] = getInstantContactTorque(b);
		if (computeContactForceField) { contactForceField[id] = getInstantContactForceField(b); }
	}
	tAccu = 0;
}

void TimeAverager::action()
{
	const Real& dt = scene->dt;
	scene->forces.sync();
	FOREACH(Body::id_t id, ids)
	{
		const shared_ptr<Body>& b = Body::byId(id, scene);
		if (!b) continue;
		pos[id]           = updateAverage(pos[id], b->state->pos, dt);
		vel[id]           = updateAverage(vel[id], b->state->vel, dt);
		angVel[id]        = updateAverage(angVel[id], b->state->angVel, dt);
		force[id]         = updateAverage(force[id], scene->forces.getForce(id), dt);
		torque[id]        = updateAverage(torque[id], scene->forces.getTorque(id), dt);
		nbContact[id]     = updateAverage(nbContact[id], getInstantNbContact(b), dt);
		contactForce[id]  = updateAverage(contactForce[id], getInstantContactForce(b), dt);
		contactTorque[id] = updateAverage(contactTorque[id], getInstantContactTorque(b), dt);
		if (computeContactForceField) {
			// Can't average directly vector<Vector3r>, so we do it element by element
			vector<Vector3r> instantContactForceField = getInstantContactForceField(b);
			for (size_t i = 0; i < contactForceField[id].size(); ++i) {
				contactForceField[id][i] = updateAverage(contactForceField[id][i], instantContactForceField[i], dt);
			}
		}
	}
	tAccu += dt;
}


// PRIVATE METHODS

template <typename T> T TimeAverager::updateAverage(const T& averagedVal, const T& instantVal, const Real& dt) const
{
	if (tAccu + dt == 0) return averagedVal;
	return (tAccu * averagedVal + dt * instantVal) / (tAccu + dt);
}

template <typename T> T TimeAverager::getValueFromMap(const boost::unordered_map<Body::id_t, T>& mapObject, Body::id_t id) const
{
	auto it = mapObject.find(id);
	if (it == mapObject.end()) { throw std::runtime_error("Particle ID not found in map"); }
	return it->second;
}

Real TimeAverager::getInstantNbContact(const shared_ptr<Body>& b) const
{
	// Get number of contacts from size of interaction map object and convert to floating point
	return b->intrs.size();
}

Vector3r TimeAverager::getInstantContactForce(const shared_ptr<Body>& b) const
{
	Vector3r instantContactForce = Vector3r::Zero();
	for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
		const shared_ptr<Interaction>& I = (*it).second;
		if (!I->isReal()) continue;
		NormShearPhys* phys = YADE_CAST<NormShearPhys*>(I->phys.get());
		instantContactForce += phys->shearForce + phys->normalForce;
	}
	return instantContactForce;
}

Vector3r TimeAverager::getInstantContactTorque(const shared_ptr<Body>& b) const
{
	Vector3r instantContactTorque = Vector3r::Zero();
	for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
		const shared_ptr<Interaction>& I = (*it).second;
		if (!I->isReal()) continue;
		NormShearPhys*         phys        = YADE_CAST<NormShearPhys*>(I->phys.get());
		GenericSpheresContact* geom        = YADE_CAST<GenericSpheresContact*>(I->geom.get());
		Vector3r               relativePos = geom->contactPoint - b->state->pos;
		Vector3r               forceSum    = phys->shearForce + phys->normalForce;
		instantContactTorque += relativePos.cross(forceSum);
	}
	return instantContactTorque;
}

vector<Vector3r> TimeAverager::getInstantContactForceField(const shared_ptr<Body>& b) const
{
	vector<Vector3r> instantContactForceField(grid.size(), Vector3r::Zero());
	// Loop on interaction map
	for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
		// Retrieve Interaction object from the current map pair
		const shared_ptr<Interaction>& I = (*it).second;
		if (!I->isReal()) continue;
		// Get  physics and geom object from the interaction
		NormShearPhys*         phys = YADE_CAST<NormShearPhys*>(I->phys.get());
		GenericSpheresContact* geom = YADE_CAST<GenericSpheresContact*>(I->geom.get());
		// Spread the contact force on each point of the grid, according to the distance between the contact point
		// and the grid point. Use a gaussian kernel to smooth the field.
		for (size_t i = 0; i < grid.size(); i++) {
			Real distanceToContact = (b->state->pos + grid[i] - geom->contactPoint).norm();
			Real weight            = exp(-pow(distanceToContact, 2) / (2 * pow(sigma, 2))) / (sqrt(2 * Mathr::PI) * sigma);
			instantContactForceField[i] += (phys->shearForce + phys->normalForce) * weight;
		}
	}
	return instantContactForceField;
}

} // namespace yade