/*****************************************************************************
*  2023 © DLH van der Haven, dannyvdhaven@gmail.com, University of Cambridge *
*                                                                            *
*  For details, please see van der Haven et al., "A physically consistent    *
*   Discrete Element Method for arbitrary shapes using Volume-interacting    *
*   Level Sets", Comput. Methods Appl. Mech. Engrg., 414 (116165):1-21       *
*   https://doi.org/10.1016/j.cma.2023.116165                                *
*  This project has been financed by Novo Nordisk A/S (Bagsværd, Denmark).   *
*                                                                            *
*  This program is licensed under GNU GPLv2, see file LICENSE for details.   *
*****************************************************************************/

#include "VolumeElasticLaw2.hpp"
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/levelSet/VolumeGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((Law2_VolumeGeom_FrictPhys_Elastic)(Law2_VolumeGeom_ViscoFrictPhys_Elastic));

Real Law2_VolumeGeom_FrictPhys_Elastic::getPlasticDissipation() const { return (Real)plasticDissipation; }

void Law2_VolumeGeom_FrictPhys_Elastic::initPlasticDissipation(Real initVal)
{
	plasticDissipation.reset();
	plasticDissipation += initVal;
}

CREATE_LOGGER(Law2_VolumeGeom_FrictPhys_Elastic);
bool Law2_VolumeGeom_FrictPhys_Elastic::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	int id1 = contact->getId1(), id2 = contact->getId2();

	VolumeGeom* geom = static_cast<VolumeGeom*>(ig.get());
	FrictPhys*  phys = static_cast<FrictPhys*>(ip.get());
	if (geom->penetrationVolume < 0) {
		if (neverErase) {
			phys->shearForce  = Vector3r::Zero();
			phys->normalForce = Vector3r::Zero();
		} else
			return false;
	}
	Real& un          = geom->penetrationVolume;
	phys->normalForce = phys->kn * math::pow(math::max(un, (Real)0), volumePower) * geom->normal;

	Vector3r&       shearForce = geom->rotate(phys->shearForce);
	const Vector3r& shearDisp  = geom->shearIncrement();
	shearForce -= phys->ks * shearDisp;
	Real maxFs = phys->normalForce.squaredNorm() * math::pow(phys->tangensOfFrictionAngle, 2);

	if (!scene->trackEnergy && !traceEnergy) { //Update force but don't compute energy terms (see below))
		// PFC3d SlipModel, is using friction angle. CoulombCriterion
		if (shearForce.squaredNorm() > maxFs) {
			Real ratio = sqrt(maxFs) / shearForce.norm();
			shearForce *= ratio;
		}
	} else {
		// Almost the same with additional Vector3r instatinated for energy tracing, duplicated block
		// to make sure there is no cost for the instanciation of the vector when traceEnergy==false
		if (shearForce.squaredNorm() > maxFs) {
			Real     ratio      = sqrt(maxFs) / shearForce.norm();
			Vector3r trialForce = shearForce; // Store previous force for definition of plastic slip
			// Define the plastic work input and increment the total plastic energy dissipated
			shearForce *= ratio;
			Real dissip = ((1 / phys->ks) * (trialForce - shearForce)) /*plastic disp*/.dot(shearForce) /*active force*/;
			if (traceEnergy) plasticDissipation += dissip;
			else if (dissip > 0)
				scene->energy->add(dissip, "plastDissip", plastDissipIx, /*reset*/ false);
		}
	}
	if (!scene->isPeriodic) { // For non-periodic simulations only
		State* de1 = Body::byId(id1, scene)->state.get();
		State* de2 = Body::byId(id2, scene)->state.get();
		applyForceAtContactPoint(-phys->normalForce - shearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position);
	} else { // The general case
		Vector3r shift2 = scene->cell->hSize * contact->cellDist.cast<Real>();
		State*   de1    = Body::byId(id1, scene)->state.get();
		State*   de2    = Body::byId(id2, scene)->state.get();
		applyForceAtContactPoint(-phys->normalForce - shearForce, geom->contactPoint, id1, de1->se3.position, id2, de2->se3.position + shift2);
	}
	return true;
}

bool Law2_VolumeGeom_ViscoFrictPhys_Elastic::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	VolumeGeom*     geom = static_cast<VolumeGeom*>(ig.get());
	ViscoFrictPhys* phys = static_cast<ViscoFrictPhys*>(ip.get());
	if (shearCreep) {
		const Real& dt = scene->dt;
		geom->rotate(phys->creepedShear);
		phys->creepedShear += creepStiffness * phys->ks * (phys->shearForce - phys->creepedShear) * dt / viscosity;
		phys->shearForce -= phys->ks * ((phys->shearForce - phys->creepedShear) * dt / viscosity);
	}
	return Law2_VolumeGeom_FrictPhys_Elastic::go(ig, ip, contact);
}

} // namespace yade
