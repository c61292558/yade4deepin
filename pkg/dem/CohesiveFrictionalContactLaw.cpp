/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@imag.fr>         *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "CohesiveFrictionalContactLaw.hpp"
#include <lib/base/LoggingUtils.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((CohesiveFrictionalContactLaw)(Law2_ScGeom6D_CohFrictPhys_CohesionMoment)(CohFrictMat)(CohFrictPhys)(Ip2_CohFrictMat_CohFrictMat_CohFrictPhys));
CREATE_LOGGER(Law2_ScGeom6D_CohFrictPhys_CohesionMoment);
CREATE_LOGGER(Ip2_CohFrictMat_CohFrictMat_CohFrictPhys);

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::getPlasticDissipation() const { return (Real)plasticDissipation; }
void Law2_ScGeom6D_CohFrictPhys_CohesionMoment::initPlasticDissipation(Real initVal)
{
	plasticDissipation.reset();
	plasticDissipation += initVal;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::normElastEnergy()
{
	Real normEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { normEnergy += 0.5 * (phys->normalForce.squaredNorm() / phys->kn); }
	}
	return normEnergy;
}
Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::shearElastEnergy()
{
	Real shearEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { shearEnergy += 0.5 * (phys->shearForce.squaredNorm() / phys->ks); }
	}
	return shearEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::bendingElastEnergy()
{
	Real bendingEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { bendingEnergy += 0.5 * (phys->moment_bending.squaredNorm() / phys->kr); }
	}
	return bendingEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::twistElastEnergy()
{
	Real twistEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { twistEnergy += 0.5 * (phys->moment_twist.squaredNorm() / phys->ktw); }
	}
	return twistEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::totalElastEnergy()
{
	Real totalEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) {
			totalEnergy += 0.5 * (phys->normalForce.squaredNorm() / phys->kn);
			totalEnergy += 0.5 * (phys->shearForce.squaredNorm() / phys->ks);
			totalEnergy += 0.5 * (phys->moment_bending.squaredNorm() / phys->kr);
			totalEnergy += 0.5 * (phys->moment_twist.squaredNorm() / phys->ktw);
		}
	}
	return totalEnergy;
}


void CohesiveFrictionalContactLaw::action()
{
	if (!functor) functor = shared_ptr<Law2_ScGeom6D_CohFrictPhys_CohesionMoment>(new Law2_ScGeom6D_CohFrictPhys_CohesionMoment);
	functor->always_use_moment_law = always_use_moment_law;
	functor->shear_creep           = shear_creep;
	functor->twist_creep           = twist_creep;
	functor->creep_viscosity       = creep_viscosity;
	functor->scene                 = scene;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		functor->go(I->geom, I->phys, I.get());
	}
}


bool Law2_ScGeom6D_CohFrictPhys_CohesionMoment::checkPlasticity(ScGeom6D* geom, CohFrictPhys* phys, Real& Fn, bool computeMoment)
{
	// if one force/torque breaks, skip computing the others since the contact properties have changed.
	// Update thresholds and compute again, with some recursive calls in section 1. Max recursion depth is 1.

	//  1. ____________________ Compare force components to their max value  _________________
	if (phys->cohesionBroken and Fn < 0) return false; // lost contact

	const bool brittle = phys->fragile and not phys->cohesionBroken;
	Real       Fs = 0, maxFs = 0, scalarRoll = 0, maxRoll = 0, scalarTwist = 0, maxTwist = 0;

	// check tensile force magnitude and max plastic disp for normal direction
	bool breaks = (-Fn) > phys->normalAdhesion or (phys->unpMax >= 0 and -(geom->penetrationDepth + phys->normalAdhesion / phys->kn) > phys->unpMax);
	// no need to compute the rest if the contact breaks
	if ((brittle or phys->unpMax >= 0) and breaks) {
		phys->SetBreakingState(always_use_moment_law);
		if (geom->penetrationDepth < 0)
			return false; // FIXME: for fragile behavior the dissipated energy is not correct (possibly), it should include the elastic energy lost in breakage
		else
			return checkPlasticity(geom, phys, Fn, computeMoment);
	}
	// check max shear force
	if (phys->ks > 0) {
		maxFs = phys->shearAdhesion;
		if ((not phys->cohesionDisablesFriction) || maxFs == 0) maxFs += Fn * phys->tangensOfFrictionAngle;
		maxFs = math::max((Real)0, maxFs);
		// compare squared norm to avoid sqrt
		Fs     = phys->shearForce.squaredNorm();
		maxFs  = maxFs * maxFs;
		breaks = Fs > maxFs;
		if (brittle and breaks) {
			phys->SetBreakingState(always_use_moment_law);
			if (geom->penetrationDepth < 0) return false;
			else
				return checkPlasticity(geom, phys, Fn, computeMoment);
		}
	}
	// check torques
	if (computeMoment) {
		// rolling torque
		if (phys->kr > 0 and (phys->maxRollPl >= 0. or phys->rollingAdhesion >= 0)) {
			maxRoll = phys->rollingAdhesion;
			if ((not phys->cohesionDisablesFriction) || maxRoll == 0) maxRoll += phys->maxRollPl * Fn;
			maxRoll    = math::max((Real)0, maxRoll);
			maxRoll    = maxRoll * maxRoll;
			scalarRoll = phys->moment_bending.squaredNorm();
			breaks     = scalarRoll > maxRoll and phys->rollingAdhesion > 0;
			if (brittle and breaks) {
				phys->SetBreakingState(always_use_moment_law);
				if (geom->penetrationDepth < 0) return false;
				else
					return checkPlasticity(geom, phys, Fn, computeMoment);
			}
		}
		// twisting torque
		if (phys->ktw > 0 and (phys->maxTwistPl >= 0. or phys->twistingAdhesion >= 0)) {
			maxTwist = phys->twistingAdhesion;
			if (maxTwist == 0 or not phys->cohesionDisablesFriction) maxTwist += phys->maxTwistPl * Fn;
			maxTwist    = math::max((Real)0, maxTwist);
			maxTwist    = maxTwist * maxTwist;
			scalarTwist = phys->moment_twist.squaredNorm();
			breaks      = scalarTwist > maxTwist and phys->twistingAdhesion > 0;
			if (brittle and breaks) {
				phys->SetBreakingState(always_use_moment_law);
				if (geom->penetrationDepth < 0) return false;
				else
					return checkPlasticity(geom, phys, Fn, computeMoment);
			}
		}
	}
	//  2. ____________________ Correct the forces and increment dissipated energy  _________________
	if ((-Fn) > phys->normalAdhesion) { //normal plasticity
		if (scene->trackEnergy || traceEnergy) {
			Real dissipated = ((1 / phys->kn) * (-Fn - phys->normalAdhesion)) * phys->normalAdhesion;
			plasticDissipation += dissipated;
			if (scene->trackEnergy) scene->energy->add(dissipated, "normalDissip", normalDissipIx, /*reset*/ false);
		}
		Fn        = -phys->normalAdhesion;
		phys->unp = geom->penetrationDepth + phys->normalAdhesion / phys->kn;
	}
	phys->normalForce = Fn * geom->normal;

	if (Fs > maxFs) { //Plasticity condition on shear force
		Real     ratio      = math::sqrt(maxFs / Fs);
		Vector3r trialForce = phys->shearForce;
		phys->shearForce *= ratio;
		if (scene->trackEnergy || traceEnergy) {
			Real sheardissip = ((1 / phys->ks) * (trialForce - phys->shearForce)) /*plastic disp*/.dot(phys->shearForce) /*active force*/;
			if (sheardissip > 0) {
				plasticDissipation += sheardissip;
				if (scene->trackEnergy) scene->energy->add(sheardissip, "shearDissip", shearDissipIx, /*reset*/ false);
			}
		}
	}
	if (computeMoment) {
		// limit rolling moment to the plastic value, if required
		if (scalarRoll > maxRoll) { // fix maximum rolling moment
			Real ratio = math::sqrt(maxRoll / scalarRoll);
			phys->moment_bending *= ratio;
			if (scene->trackEnergy) {
				Real bendingdissip = ((1 / phys->kr) * (scalarRoll - maxRoll) * maxRoll) /*active force*/;
				if (bendingdissip > 0) scene->energy->add(bendingdissip, "bendingDissip", bendingDissipIx, /*reset*/ false);
			}
		}
		// limit twisting moment to the plastic value, if required
		if (scalarTwist > maxTwist) { // fix maximum rolling moment
			Real ratio = math::sqrt(maxTwist / scalarTwist);
			phys->moment_twist *= ratio;
			if (scene->trackEnergy) {
				Real twistdissip = ((1 / phys->ktw) * (scalarTwist - maxTwist) * maxTwist) /*active force*/;
				if (twistdissip > 0) scene->energy->add(twistdissip, "twistDissip", twistDissipIx, /*reset*/ false);
			}
		}
	}
	return true;
}


void Law2_ScGeom6D_CohFrictPhys_CohesionMoment::setElasticForces(
        shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I, bool computeMoment, Real& Fn, const Real& dt)
{
	State*        de1  = Body::byId(I->getId1(), scene)->state.get();
	State*        de2  = Body::byId(I->getId2(), scene)->state.get();
	ScGeom6D*     geom = YADE_CAST<ScGeom6D*>(ig.get());
	CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(ip.get());
	// NORMAL
	Fn = phys->kn * (geom->penetrationDepth - phys->unp);

	// SHEAR
	// 	// update orientation
	// 	Vector3r&       shearForce = geom->rotate(phys->shearForce);
	// increment
	const Vector3r& dus = geom->shearIncrement();
	phys->shearForce -= phys->ks * dus;

	// TORQUES
	if (computeMoment) {
		if (!useIncrementalForm) {
			if (phys->ktw > 0) phys->moment_twist = (geom->getTwist() * phys->ktw) * geom->normal;
			else
				phys->moment_twist = Vector3r::Zero();
			if (phys->kr > 0) phys->moment_bending = geom->getBending() * phys->kr;
			else
				phys->moment_bending = Vector3r::Zero();
		} else { // Use incremental formulation to compute moment_twis and moment_bending (no twist_creep is applied)
			Vector3r relAngVel = geom->getRelAngVel(de1, de2, dt);
			// *** Bending ***//
			if (phys->kr > 0) {
				Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel) * geom->normal; // keep only the bending part
				// incremental formulation for the bending moment (as for the shear part)
				phys->moment_bending = phys->moment_bending - phys->kr * relAngVelBend * dt;
			} else
				phys->moment_bending = Vector3r::Zero();
			// ----------------------------------------------------------------------------------------
			// *** Torsion ***//
			if (phys->ktw > 0) {
				Vector3r relAngVelTwist = geom->normal.dot(relAngVel) * geom->normal;
				Vector3r relRotTwist    = relAngVelTwist * dt; // component of relative rotation along n
				phys->moment_twist      = phys->moment_twist - phys->ktw * relRotTwist;
			} else
				phys->moment_twist = Vector3r::Zero();
		}
	}
}

void Law2_ScGeom6D_CohFrictPhys_CohesionMoment::checkConsistency(const shared_ptr<CohFrictPhys> phys, Body::id_t id1, Body::id_t id2) const
{
	if (phys->shearAdhesion > 0 and (phys->maxRollPl < 0. or phys->rollingAdhesion <= 0))
		LOG_WARN_ONCE(
		        "the max shear force includes some adhesion/cohesion but the rolling moment does not (because it is either purely elastic or purely "
		        "frictional) for interaction ("
		        << id1 << "-" << id2 << "), if it is not intentional you may want to check input parameters.");
	if ((phys->maxRollPl >= 0. or phys->rollingAdhesion > 0) and (phys->maxTwistPl < 0. and phys->ktw > 0))
		LOG_WARN_ONCE(
		        "rolling moment is bounded but twisting moment is not (interaction "
		        << id1 << "-" << id2 << "), check etaRoll and etaTwist of the materials if it is not intentional.");
	if (phys->maxRollPl < 0. and (phys->maxTwistPl >= 0. or phys->twistingAdhesion > 0))
		LOG_WARN_ONCE(
		        "twisting moment is bounded but rolling moment is not (interaction "
		        << id1 << "-" << id2 << "), check etaRoll and etaTwist of the materials if it is not intentional.");
	if ((phys->maxRollPl >= 0. or phys->rollingAdhesion >= 0 or phys->maxTwistPl >= 0. or phys->twistingAdhesion >= 0) and not useIncrementalForm)
		LOG_WARN_ONCE("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not be applied "
		              "correctly to torques (the total formulation will not reproduce irreversibility).");
	if ((phys->maxRollPl != 0 or phys->maxTwistPl != 0 or phys->rollingAdhesion != 0 or phys->twistingAdhesion != 0) and not phys->momentRotationLaw)
		LOG_WARN_ONCE(
		        "Interaction " << id1 << "-" << id2
		                       << " has some parameters related to contact moment defined but it will be ignored because i.phys.momentRotationLaw is "
		                          "False. Make sure momentRotationLaw is True in the material class if moments are needed.");
}

bool Law2_ScGeom6D_CohFrictPhys_CohesionMoment::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real&   dt              = scene->dt;
	const int&    id1             = contact->getId1();
	const int&    id2             = contact->getId2();
	State*        de1             = Body::byId(id1, scene)->state.get();
	State*        de2             = Body::byId(id2, scene)->state.get();
	ScGeom6D*     geom            = YADE_CAST<ScGeom6D*>(ig.get());
	CohFrictPhys* phys            = YADE_CAST<CohFrictPhys*>(ip.get());
	Vector3r&     shearForceFirst = phys->shearForce;
	bool          computeMoment   = phys->momentRotationLaw and ((!phys->cohesionBroken) or always_use_moment_law);

	// FIXME: move this to Ip2::go and put the code block in a separate function)
	if (consistencyCheck and phys->momentRotationLaw) {
		checkConsistency(YADE_PTR_CAST<CohFrictPhys>(ip), id1, id2);
		consistencyCheck = false;
	}

	if (contact->isFresh(scene)) shearForceFirst = Vector3r::Zero();

	//  ____________________ Creep current forces if necessary _________________
	if (shear_creep) shearForceFirst -= phys->ks * (shearForceFirst * dt / creep_viscosity);
	if (twist_creep) {
		if (!useIncrementalForm) {
			Real        viscosity_twist     = creep_viscosity * math::pow((2 * math::min(geom->radius1, geom->radius2)), 2) / 16.0;
			Real        angle_twist_creeped = geom->getTwist() * (1 - dt / viscosity_twist);
			Quaternionr q_twist(AngleAxisr(geom->getTwist(), geom->normal));
			Quaternionr q_twist_creeped(AngleAxisr(angle_twist_creeped, geom->normal));
			Quaternionr q_twist_delta(q_twist_creeped * q_twist.conjugate());
			geom->twistCreep = geom->twistCreep * q_twist_delta;
		} else
			LOG_WARN_ONCE("Law2_ScGeom6D_CohFrictPhys_CohesionMoment: no twis creep is included if the incremental formulation is used.");
	}

	//  ____________________  Update orientation ___________________________________
	geom->rotate(phys->shearForce);
	if (computeMoment and useIncrementalForm) {
		if (phys->kr > 0) geom->rotate(phys->moment_bending);
		if (phys->ktw > 0) geom->rotate(phys->moment_twist);
	}

	//  ____________________ Linear elasticity giving "trial" force/torques _________________

	Real Fn = 0;
	setElasticForces(ig, ip, contact, computeMoment, Fn, dt);

	//  ____________________ Failure and plasticity _________________________________________

	bool alive = checkPlasticity(geom, phys, Fn, computeMoment);

	// ____________________  Delete interaction if no contact  ______________________________
	if (not alive) { return false; }
	//  ____________________  Apply the force _______________________________________________
	applyForceAtContactPoint(
	        -phys->normalForce - phys->shearForce,
	        geom->contactPoint,
	        id1,
	        de1->se3.position,
	        id2,
	        de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist) : Vector3r::Zero()));

	if (computeMoment) {
		// Apply moments now
		Vector3r moment = phys->moment_twist + phys->moment_bending;
		scene->forces.addTorque(id1, -moment);
		scene->forces.addTorque(id2, moment);
	}
	return true;
}

bool Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, CohFrictPhys* contactPhysics)
{
	if (not contactPhysics) {
		contactPhysics = YADE_DYN_CAST<CohFrictPhys*>(interaction->phys.get());
		if (not contactPhysics) {
			LOG_WARN("Invalid type of interaction, cohesion not set");
			return false;
		}
	}
	// if breaks
	if ((not cohesive) and not contactPhysics->cohesionBroken) {
		contactPhysics->SetBreakingState(true);
		return true;
	}
	// else bond
	if (not scene) scene = Omega::instance().getScene().get();
	auto       b1   = Body::byId(interaction->getId1(), scene);
	auto       b2   = Body::byId(interaction->getId2(), scene);
	auto       mat1 = static_cast<CohFrictMat*>(b1->material.get());
	auto       mat2 = static_cast<CohFrictMat*>(b2->material.get());
	const auto geom = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	// determine adhesion terms with matchmakers or default formulae
	const auto normalCohPreCalculated   = (normalCohesion) ? (*normalCohesion)(b1->id, b2->id) : math::min(mat1->normalCohesion, mat2->normalCohesion);
	const auto shearCohPreCalculated    = (shearCohesion) ? (*shearCohesion)(b1->id, b2->id) : math::min(mat1->shearCohesion, mat2->shearCohesion);
	const auto rollingCohPreCalculated  = (rollingCohesion) ? (*rollingCohesion)(b1->id, b2->id) : normalCohPreCalculated;
	const auto twistingCohPreCalculated = (twistingCohesion) ? (*twistingCohesion)(b1->id, b2->id) : shearCohPreCalculated;
	// assign adhesions
	contactPhysics->cohesionBroken = false;
	contactPhysics->normalAdhesion = normalCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
	contactPhysics->shearAdhesion  = shearCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
	if (contactPhysics->momentRotationLaw) {
		// the max stress in pure bending is 4*M/πr^3 = 4*M/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/4)*normalAdhesion
		contactPhysics->rollingAdhesion = 0.25 * rollingCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
		// the max shear stress in pure twisting is 2*Mt/πr^3 = 2*Mt/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/2)*shearAdhesion
		contactPhysics->twistingAdhesion = 0.5 * twistingCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
	}
	geom->initRotations(*(b1->state), *(b2->state));
	contactPhysics->fragile      = (mat1->fragile || mat2->fragile);
	contactPhysics->initCohesion = false;
	return true;
}


void Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::go(
        const shared_ptr<Material>& b1 // CohFrictMat
        ,
        const shared_ptr<Material>& b2 // CohFrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	CohFrictMat* mat1 = static_cast<CohFrictMat*>(b1.get());
	CohFrictMat* mat2 = static_cast<CohFrictMat*>(b2.get());
	ScGeom6D*    geom = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	//Create cohesive interractions only once
	if (setCohesionNow && cohesionDefinitionIteration == -1) cohesionDefinitionIteration = scene->iter;
	if (setCohesionNow && cohesionDefinitionIteration != -1 && cohesionDefinitionIteration != scene->iter) {
		cohesionDefinitionIteration = -1;
		setCohesionNow              = 0;
	}

	if (geom) {
		if (!interaction->phys) {
			interaction->phys            = shared_ptr<CohFrictPhys>(new CohFrictPhys());
			CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
			Real          Ea             = mat1->young;
			Real          Eb             = mat2->young;
			Real          Va             = mat1->poisson;
			Real          Vb             = mat2->poisson;
			Real          Da             = geom->radius1;
			Real          Db             = geom->radius2;
			Real          fa             = mat1->frictionAngle;
			Real          fb             = mat2->frictionAngle;
			Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			Real          frictionAngle  = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(mat1->id, mat2->id, fa, fb);

			// harmonic average of alphas parameters
			Real AlphaKr, AlphaKtw;
			if (mat1->alphaKr && mat2->alphaKr) AlphaKr = 2.0 * mat1->alphaKr * mat2->alphaKr / (mat1->alphaKr + mat2->alphaKr);
			else
				AlphaKr = 0;
			if (mat1->alphaKtw && mat2->alphaKtw) AlphaKtw = 2.0 * mat1->alphaKtw * mat2->alphaKtw / (mat1->alphaKtw + mat2->alphaKtw);
			else
				AlphaKtw = 0;

			Real Ks;
			if (Va && Vb)
				Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb
				        / (Ea * Da * Va + Eb * Db * Vb); //harmonic average of two stiffnesses with ks=V*kn for each sphere
			else
				Ks = 0;

			contactPhysics->kn                     = Kn;
			contactPhysics->ks                     = Ks;
			contactPhysics->kr                     = Da * Db * Ks * AlphaKr;
			contactPhysics->ktw                    = Da * Db * Ks * AlphaKtw;
			contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
			contactPhysics->momentRotationLaw      = (mat1->momentRotationLaw && mat2->momentRotationLaw);
			if (contactPhysics->momentRotationLaw) {
				contactPhysics->maxRollPl  = min(mat1->etaRoll * Da, mat2->etaRoll * Db);
				contactPhysics->maxTwistPl = min(mat1->etaTwist * Da, mat2->etaTwist * Db);
			} else {
				contactPhysics->maxRollPl = contactPhysics->maxTwistPl = 0;
			}
			if ((setCohesionOnNewContacts || setCohesionNow) && mat1->isCohesive && mat2->isCohesive) {
				setCohesion(interaction, true, contactPhysics);
			}
		} else { // !isNew, but if setCohesionNow, all contacts are initialized like if they were newly created
			CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
			if ((setCohesionNow or contactPhysics->initCohesion) and (mat1->isCohesive && mat2->isCohesive))
				setCohesion(interaction, true, contactPhysics);
		}
	}
};

} // namespace yade
