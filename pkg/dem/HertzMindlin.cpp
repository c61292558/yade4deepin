// 2010 © Chiara Modenese <c.modenese@gmail.com>

#include "HertzMindlin.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((MindlinPhys)(Ip2_FrictMat_FrictMat_MindlinPhys)(Law2_ScGeom_MindlinPhys_MindlinDeresiewitz)(Law2_ScGeom_MindlinPhys_HertzWithLinearShear)(
        Law2_ScGeom_MindlinPhys_Mindlin)(MindlinCapillaryPhys)(Ip2_FrictMat_FrictMat_MindlinCapillaryPhys));

Real Law2_ScGeom_MindlinPhys_Mindlin::getfrictionDissipation() const { return (Real)frictionDissipation; }
Real Law2_ScGeom_MindlinPhys_Mindlin::getshearEnergy() const { return (Real)shearEnergy; }
Real Law2_ScGeom_MindlinPhys_Mindlin::getnormDampDissip() const { return (Real)normDampDissip; }
Real Law2_ScGeom_MindlinPhys_Mindlin::getshearDampDissip() const { return (Real)shearDampDissip; }


/* Functions to calculate velocity-dependent coefficient of restitution as in [Brilliantov1996]_ and [Mueller2011]_ */

// Function to calculate the restitution coefficient as a function of normalised velocity for the viscous damping model of [Brilliantov1996]_ using Pade approximation, as in [Mueller2011]_
Real restitutionCoefficient(const Real v_star)
{
	// 3/6 coefficients
	const Real a_i[] = { 1.0, 1.07232, 0.574198, 0.141552 };
	const Real b_i[] = { 1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585 };

	// Initialize sum to 0
	Real A = 0.0, B = 0.0, n = 0.0;

	for (auto& i : a_i) {
		A += i * math::pow(v_star, n);
		n += 1.0;
	}

	n = 0.0;
	for (auto& i : b_i) {
		B += i * math::pow(v_star, n);
		n += 1.0;
	}

	return A / B;
}

// Function to calculate derivative used in Newton-Raphson iterations
Real restitutionCoefficientDeriv(const Real v_star)
{
	// 3/6 coefficients
	const Real a_i[] = { 1.0, 1.07232, 0.574198, 0.141552 };
	const Real b_i[] = { 1.0, 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585 };

	// Initialize sums to 0
	Real A = 0.0, B = 0.0, dA = 0.0, dB = 0.0, n = 0.0;

	for (auto& i : a_i) {
		A += i * math::pow(v_star, n);
		n += 1.0;
	}

	n = 0.0;
	for (auto& i : b_i) {
		B += i * math::pow(v_star, n);
		n += 1.0;
	}

	// 3/6 coefficients Derivative
	const Real da_i[] = { 1.07232, 0.574198, 0.141552 };
	const Real db_i[] = { 1.07232, 1.72765, 1.37842, 1.19449, 0.467273, 0.235585 };

	n = 0.0;
	for (auto& i : da_i) {
		dA += (n + 1.0) * i * math::pow(v_star, n);
		n += 1.0;
	}

	n = 0.0;
	for (auto& i : db_i) {
		dB += (n + 1.0) * i * math::pow(v_star, n);
		n += 1.0;
	}

	return (dA * B - A * dB) / (B * B);
}


// Function to calculate vstar as a function of the coefficient of restitution according to [Mueller2011]_
Real getVstar(const Real en)
{
	const int max_iter = 1000; // Maximum number of iterations
	int       i        = 0;
	Real      xr       = 0.5; // Initial guess

	// Check if en is in the range [0,1]
	if (en == 1.0) return 0.0;
	else if (en > 1.0)
		throw std::runtime_error("getVstar: en > 1. Restitution coefficient must be between 0 and 1.");
	else if (en < 0.0)
		throw std::runtime_error("getVstar: en < 0. Restitution coefficient must be between 0 and 1.");

	// Newton Rapshon method
	for (i = 0; i < max_iter; ++i) {
		// Check if root was found
		if (math::abs(restitutionCoefficient(xr) - en) <= std::numeric_limits<Real>::epsilon()) break;

		// New root approximation
		xr -= (restitutionCoefficient(xr) - en) / restitutionCoefficientDeriv(xr);
	}

	// Check if maximum number of iterations was reached
	if (i >= max_iter) {
		std::cerr << "WARNING: getVstar: Maximum number of iterations reached. Root not found." << std::endl;
		std::cerr << "WARNING: getVstar: | f(x) - targuet | = " << math::fabs(restitutionCoefficient(xr) - en) << std::endl;
	}

	return xr;
}


/******************** Ip2_FrictMat_FrictMat_MindlinPhys *******/
CREATE_LOGGER(Ip2_FrictMat_FrictMat_MindlinPhys);

void Ip2_FrictMat_FrictMat_MindlinPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return; // no updates of an already existing contact necessary
	shared_ptr<MindlinPhys> contactPhysics(new MindlinPhys());
	interaction->phys = contactPhysics;
	const auto mat1   = YADE_CAST<FrictMat*>(b1.get());
	const auto mat2   = YADE_CAST<FrictMat*>(b2.get());

	/* from interaction physics */
	const Real Ea = mat1->young;
	const Real Eb = mat2->young;
	const Real Va = mat1->poisson;
	const Real Vb = mat2->poisson;
	const Real fa = mat1->frictionAngle;
	const Real fb = mat2->frictionAngle;


	/* from interaction geometry */
	const auto scg = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	const Real Da  = scg->refR1 > 0 ? scg->refR1 : scg->refR2;
	const Real Db  = scg->refR2 > 0 ? scg->refR2 : scg->refR1;
	//Vector3r normal=scg->normal;        //The variable set but not used


	/* calculate stiffness coefficients */
	const Real Ga = Ea / (2 * (1 + Va));
	const Real Gb = Eb / (2 * (1 + Vb));
	const Real G  = 1.0 / ((2 - Va) / Ga + (2 - Vb) / Gb); // effective shear modulus
	//	const Real V             = (Va + Vb) / 2;                                                           // average of poisson's ratio
	const Real E             = Ea * Eb / ((1. - math::pow(Va, 2)) * Eb + (1. - math::pow(Vb, 2)) * Ea); // effective Young modulus
	const Real R             = Da * Db / (Da + Db);                                                     // equivalent radius
	const Real Rmean         = (Da + Db) / 2.;                                                          // mean radius
	const Real Kno           = 4. / 3. * E * sqrt(R);                                                   // coefficient for normal stiffness
	const Real Kso           = 8 * sqrt(R) * G;                                                         // coefficient for shear stiffness
	const Real frictionAngle = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(mat1->id, mat2->id, mat1->frictionAngle, mat2->frictionAngle);

	const Real Adhesion = 4. * Mathr::PI * R * gamma; // calculate adhesion force as predicted by DMT theory

	/* pass values calculated from above to MindlinPhys */
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	//contactPhysics->prevNormal = scg->normal; // used to compute relative rotation
	contactPhysics->kno           = Kno; // this is just a coeff
	contactPhysics->kso           = Kso; // this is just a coeff
	contactPhysics->adhesionForce = Adhesion;

	contactPhysics->kr        = krot;
	contactPhysics->ktw       = ktwist;
	contactPhysics->maxBendPl = eta * Rmean; // does this make sense? why do we take Rmean?

	/* compute viscous coefficients */
	if (en && betan) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinPhys: only one of en, betan can be specified.");
	if (es && betas) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinPhys: only one of es, betas can be specified.");

	if (vn && betan) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinPhys: only one of vn, betan can be specified.");
	if (vn && betas) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinPhys: only one of vn, betas can be specified.");

	// en or es specified
	if (en && vn) { // Velocity-dependent coefficient of restitution
		Real vstar           = getVstar((*en)(mat1->id, mat2->id));
		contactPhysics->beta = vstar * vstar * math::pow((*vn)(mat1->id, mat2->id), -0.2); // ^-1/5 = ^-2/10
		if (es) { std::cout << "Since vn is defined, the shear coefficient of restitution es will not be used." << endl; }
	} else if (en || es) {              // Constant coefficient of restitution
		const Real h1  = -6.918798; // Fitting coefficients h_i from  Table 2 - Thornton et al. (2013).
		const Real h2  = -16.41105;
		const Real h3  = 146.8049;
		const Real h4  = -796.4559;
		const Real h5  = 2928.711;
		const Real h6  = -7206.864;
		const Real h7  = 11494.29;
		const Real h8  = -11342.18;
		const Real h9  = 6276.757;
		const Real h10 = -1489.915;

		// Consider same coefficient of restitution if only one is given (en or es)
		if (!en) { en = es; }
		if (!es) { es = en; }

		const Real En     = (*en)(mat1->id, mat2->id);
		const Real Es     = (*es)(mat1->id, mat2->id);
		const Real alphan = En
		        * (h1
		           + En * (h2 + En * (h3 + En * (h4 + En * (h5 + En * (h6 + En * (h7 + En * (h8 + En * (h9 + En * h10))))))))); // Eq. (B7) from Thornton et al. (2013)
		contactPhysics->betan = (En == 1.0) ? 0
		                                    : sqrt(1.0 / (1.0 - (math::pow(1.0 + En, 2)) * exp(alphan))
		                                           - 1.0); // Eq. (B6) from Thornton et al. (2013) - This is noted as 'gamma' in their paper

		// although Thornton (2015) considered betan=betas, here we use his formulae (B6) and (B7) allowing for betas to take a different value, based on the input es
		const Real alphas     = Es * (h1 + Es * (h2 + Es * (h3 + Es * (h4 + Es * (h5 + Es * (h6 + Es * (h7 + Es * (h8 + Es * (h9 + Es * h10)))))))));
		contactPhysics->betas = (Es == 1.0) ? 0 : sqrt(1.0 / (1.0 - (math::pow(1.0 + Es, 2)) * exp(alphas)) - 1.0);

		// betan/betas specified, use that value directly
	} else { // Constant coefficient of restitution
		contactPhysics->betan = betan ? (*betan)(mat1->id, mat2->id) : 0;
		contactPhysics->betas = betas ? (*betas)(mat1->id, mat2->id) : contactPhysics->betan;
	}
}

/* Function to count the number of adhesive contacts in the simulation at each time step */
Real Law2_ScGeom_MindlinPhys_Mindlin::contactsAdhesive() // It is returning something rather than zero only if includeAdhesion is set to true
{
	Real contactsAdhesive = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		MindlinPhys* phys = dynamic_cast<MindlinPhys*>(I->phys.get());
		if (phys->isAdhesive) { contactsAdhesive += 1; }
	}
	return contactsAdhesive;
}

/* Function which returns the ratio between the number of sliding contacts to the total number at a given time */
Real Law2_ScGeom_MindlinPhys_Mindlin::ratioSlidingContacts()
{
	Real ratio(0);
	int  count(0);
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		MindlinPhys* phys = dynamic_cast<MindlinPhys*>(I->phys.get());
		if (phys->isSliding) { ratio += 1; }
		count++;
	}
	ratio /= count;
	return ratio;
}

/* Function to get the NORMAL elastic potential energy of the system */
Real Law2_ScGeom_MindlinPhys_Mindlin::normElastEnergy()
{
	Real normEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		ScGeom*      scg  = dynamic_cast<ScGeom*>(I->geom.get());
		MindlinPhys* phys = dynamic_cast<MindlinPhys*>(I->phys.get());
		if (phys) {
			if (includeAdhesion) {
				normEnergy += (math::pow(scg->penetrationDepth, 5. / 2.) * 2. / 5. * phys->kno - phys->adhesionForce * scg->penetrationDepth);
			} else {
				normEnergy += math::pow(scg->penetrationDepth, 5. / 2.) * 2. / 5. * phys->kno;
			} // work done in the normal direction. NOTE: this is the integral
		}
	}
	return normEnergy;
}

/* Function to get the adhesion energy of the system */
Real Law2_ScGeom_MindlinPhys_Mindlin::adhesionEnergy()
{
	Real adhesionEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		ScGeom*      scg  = dynamic_cast<ScGeom*>(I->geom.get());
		MindlinPhys* phys = dynamic_cast<MindlinPhys*>(I->phys.get());
		if (phys && includeAdhesion) {
			Real R       = scg->radius1 * scg->radius2 / (scg->radius1 + scg->radius2);
			Real gammapi = phys->adhesionForce / (4. * R);
			adhesionEnergy += gammapi * pow(phys->radius, 2);
		} // note that contact radius is calculated if we calculate energy components
	}
	return adhesionEnergy;
}

bool Law2_ScGeom_MindlinPhys_MindlinDeresiewitz::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	Body::id_t   id1(contact->getId1()), id2(contact->getId2());
	ScGeom*      geom = static_cast<ScGeom*>(ig.get());
	MindlinPhys* phys = static_cast<MindlinPhys*>(ip.get());
	const Real   uN   = geom->penetrationDepth;
	if (uN < 0) {
		if (neverErase) {
			phys->shearForce = phys->normalForce = Vector3r::Zero();
			phys->kn = phys->ks = 0;
			return true;
		} else {
			return false;
		}
	}
	// normal force
	Real Fn           = phys->kno * pow(uN, 3 / 2.);
	phys->normalForce = Fn * geom->normal;
	// exactly zero would not work with the shear formulation, and would give zero shear force anyway
	if (Fn == 0) return true;
	//phys->kn=3./2.*phys->kno*math::pow(uN,0.5); // update stiffness, not needed

	// contact radius
	Real R       = geom->radius1 * geom->radius2 / (geom->radius1 + geom->radius2);
	phys->radius = pow(Fn * pow(R, 3 / 2.) / phys->kno, 1 / 3.);

	// shear force: transform, but keep the old value for now
	geom->rotate(phys->usTotal);
	//Vector3r usOld=phys->usTotal;     //The variable set but not used
	Vector3r dUs = geom->shearIncrement();
	phys->usTotal -= dUs;

#if 0
	Vector3r shearIncrement;
	shearIncrement=geom->shearIncrement();
	Fs-=ks*shearIncrement;
	// Mohr-Coulomb slip
	Real maxFs2=pow(Fn,2)*pow(phys->tangensOfFrictionAngle,2);
	if(Fs.squaredNorm()>maxFs2) Fs*=sqrt(maxFs2)/Fs.norm();
#endif
	// apply forces
	Vector3r f = -phys->normalForce - phys->shearForce;
	scene->forces.addForce(id1, f);
	scene->forces.addForce(id2, -f);
	scene->forces.addTorque(id1, (geom->radius1 - .5 * geom->penetrationDepth) * geom->normal.cross(f));
	scene->forces.addTorque(id2, (geom->radius2 - .5 * geom->penetrationDepth) * geom->normal.cross(f));
	return true;
}

bool Law2_ScGeom_MindlinPhys_HertzWithLinearShear::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	Body::id_t   id1(contact->getId1()), id2(contact->getId2());
	ScGeom*      geom = static_cast<ScGeom*>(ig.get());
	MindlinPhys* phys = static_cast<MindlinPhys*>(ip.get());
	const Real   uN   = geom->penetrationDepth;
	if (uN < 0) {
		if (neverErase) {
			phys->shearForce = phys->normalForce = Vector3r::Zero();
			phys->kn = phys->ks = 0;
			return true;
		} else
			return false;
	}
	// normal force
	Real Fn           = phys->kno * pow(uN, 3 / 2.);
	phys->normalForce = Fn * geom->normal;
	//phys->kn=3./2.*phys->kno*math::pow(uN,0.5); // update stiffness, not needed

	// shear force
	Vector3r&       Fs             = geom->rotate(phys->shearForce);
	Real            ks             = nonLin > 0 ? phys->kso * math::pow(uN, 0.5) : phys->kso;
	const Vector3r& shearIncrement = geom->shearIncrement();
	Fs -= ks * shearIncrement;
	// Mohr-Coulomb slip
	Real maxFs2 = pow(Fn, 2) * pow(phys->tangensOfFrictionAngle, 2);
	if (Fs.squaredNorm() > maxFs2) Fs *= sqrt(maxFs2) / Fs.norm();

	// apply forces
	Vector3r f = -phys->normalForce - phys->shearForce; /* should be a reference returned by geom->rotate */
	assert(phys->shearForce == Fs);
	scene->forces.addForce(id1, f);
	scene->forces.addForce(id2, -f);
	scene->forces.addTorque(id1, (geom->radius1 - .5 * geom->penetrationDepth) * geom->normal.cross(f));
	scene->forces.addTorque(id2, (geom->radius2 - .5 * geom->penetrationDepth) * geom->normal.cross(f));
	return true;
}


/******************** Law2_ScGeom_MindlinPhys_Mindlin *********/
CREATE_LOGGER(Law2_ScGeom_MindlinPhys_Mindlin);

bool Law2_ScGeom_MindlinPhys_Mindlin::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real& dt = scene->dt; // get time step

	Body::id_t id1 = contact->getId1(); // get id body 1
	Body::id_t id2 = contact->getId2(); // get id body 2

	auto* de1 = Body::byId(id1, scene)->state.get();
	auto* de2 = Body::byId(id2, scene)->state.get();

	ScGeom*      scg  = static_cast<ScGeom*>(ig.get());
	MindlinPhys* phys = static_cast<MindlinPhys*>(ip.get());

	const shared_ptr<Body>& b1 = Body::byId(id1, scene);
	const shared_ptr<Body>& b2 = Body::byId(id2, scene);

	bool useDamping = (phys->betan != 0. || phys->betas != 0. || phys->beta != 0.);

#ifdef PARTIALSAT
	if (contact->isFresh(scene)) {
		phys->initD = scg->penetrationDepth; // only useful for partialsat break criteria
	}
#endif

	Real cn = 0, cs = 0;

	/****************/
	/* NORMAL FORCE */
	/****************/

	Real uN = scg->penetrationDepth; // get overlapping
	if (uN < 0) {
		if (neverErase) {
			phys->shearForce = phys->normalForce = Vector3r::Zero();
			phys->kn = phys->ks = 0;
			return true;
		} else
			return false;
	}
	/* Hertz-Mindlin's formulation (PFC)
	Note that the normal stiffness here is a secant value (so as it is cannot be used in the GSTS)
	In the first place we get the normal force and then we store kn to be passed to the GSTS */
	Real Fn = phys->kno * math::pow(uN, 1.5); // normal Force (scalar)
	if (includeAdhesion) {
		Fn -= phys->adhesionForce;   // include adhesion force to account for the effect of Van der Waals interactions
		phys->isAdhesive = (Fn < 0); // set true the bool to count the number of adhesive contacts
	}
	phys->normalForce = Fn * scg->normal; // normal Force (vector)

	if (calcEnergy) {
		Real R = scg->radius1 * scg->radius2 / (scg->radius1 + scg->radius2);
		phys->radius
		        = pow((Fn + (includeAdhesion ? phys->adhesionForce : 0.)) * pow(R, 3 / 2.) / phys->kno,
		              1 / 3.); // attribute not used anywhere, we do not need it
	}

	/*******************************/
	/* TANGENTIAL NORMAL STIFFNESS */
	/*******************************/

	phys->kn = 3. / 2. * phys->kno * math::pow(uN, 0.5); // here we store the value of kn to compute the time step

	/******************************/
	/* TANGENTIAL SHEAR STIFFNESS */
	/******************************/

	phys->ks = phys->kso * math::pow(uN, 0.5); // get tangential stiffness (this is a tangent value, so we can pass it to the GSTS)

	/************************/
	/* DAMPING COEFFICIENTS */
	/************************/

	if (useDamping) {
		Real mbar = (!b1->isDynamic() && b2->isDynamic())
		        ? de2->mass
		        : ((!b2->isDynamic() && b1->isDynamic())
		                   ? de1->mass
		                   : (de1->mass * de2->mass
		                      / (de1->mass
		                         + de2->mass))); // get equivalent mass if both bodies are dynamic, if not set it equal to the one of the dynamic body
		if (phys->betan != 0. || phys->betas != 0.) {      // Constant coefficient of restitution (see Thornton, 2015)
			Real Cn_crit = 2. * sqrt(mbar * phys->kn); // Critical damping coefficient (normal direction)
			Real Cs_crit = 2. * sqrt(mbar * phys->ks); // Critical damping coefficient (shear direction)

			cn = Cn_crit * phys->betan; // Damping normal coefficient
			cs = Cs_crit * phys->betas; // Damping tangential coefficient
			if (phys->kn < 0 || phys->ks < 0) {
				cerr << "Negative stiffness kn=" << phys->kn << " ks=" << phys->ks << " for ##" << b1->getId() << "+" << b2->getId()
				     << ", step " << scene->iter << endl;
			}
		} else if (phys->beta != 0.) {                                                   // Velocity-dependent coefficient of restitution
			const Real A = 2.0 * phys->beta * math::pow(phys->kno / mbar, -0.4) / 3; // ^-2/5
			cn           = A * phys->kn;
			cs           = A * phys->kn;
		}
	}

	/***************/
	/* SHEAR FORCE */
	/***************/

	Vector3r& shearElastic = phys->shearElastic; // reference for shearElastic force
	// Define shifts to handle periodicity
	const Vector3r shift2   = scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist) : Vector3r::Zero();
	const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(contact->cellDist) : Vector3r::Zero();
	// 1. Rotate shear force
	shearElastic            = scg->rotate(shearElastic);
	Vector3r prev_FsElastic = shearElastic; // save shear force at previous time step
	                                        // 2. Get incident velocity, get shear and normal components
	// NOTE: below, the normal component is obtained from getIncidentVel(), OTOH, the shear component computed at next line would be wrong for sphere-facet
	// and possibly other Ig types incompatible with preventGranularRatcheting=true. We thus use the precomputed shearIncrement from the Ig2, which should be always correct.
	Vector3r incidentV = scg->getIncidentVel(de1, de2, dt, shift2, shiftVel, false);
	//     Vector3r incidentV  = geom->shearIncrement()/dt;
	Vector3r incidentVn = scg->normal.dot(incidentV) * scg->normal; // contact normal velocity
	Vector3r incidentVs = scg->shearIncrement() / dt;               // contact shear velocity
	// 3. Get shear force (incrementally)
	shearElastic = shearElastic - phys->ks * (incidentVs * dt);

	/**************************************/
	/* VISCOUS DAMPING (Normal direction) */
	/**************************************/

	// normal force must be updated here before we apply the Mohr-Coulomb criterion
	if (useDamping) { // get normal viscous component
		phys->normalViscous = cn * incidentVn;
		Vector3r normTemp   = phys->normalForce - phys->normalViscous; // temporary normal force
		// viscous force should not exceed the value of current normal force, i.e. no attraction force should be permitted if particles are non-adhesive
		// if particles are adhesive, then fixed the viscous force at maximum equal to the adhesion force
		// *** enforce normal force to zero if no adhesion is permitted ***
		if (phys->adhesionForce == 0.0 || !includeAdhesion) {
			if (normTemp.dot(scg->normal) < 0.0) {
				phys->normalForce   = Vector3r::Zero();
				phys->normalViscous = phys->normalViscous
				        + normTemp; // normal viscous force is such that the total applied force is null - it is necessary to compute energy correctly!
			} else {
				phys->normalForce -= phys->normalViscous;
			}
		} else if (includeAdhesion && phys->adhesionForce != 0.0) {
			// *** limit viscous component to the max adhesive force ***
			if (normTemp.dot(scg->normal) < 0.0 && (phys->normalViscous.norm() > phys->adhesionForce)) {
				Real     normVisc       = phys->normalViscous.norm();
				Vector3r normViscVector = phys->normalViscous / normVisc;
				phys->normalViscous     = phys->adhesionForce * normViscVector;
				phys->normalForce -= phys->normalViscous;
			}
			// *** apply viscous component - in the presence of adhesion ***
			else {
				phys->normalForce -= phys->normalViscous;
			}
		}
		if (calcEnergy) { normDampDissip += phys->normalViscous.dot(incidentVn * dt); } // calc dissipation of energy due to normal damping
	}


	/*************************************/
	/* SHEAR DISPLACEMENT (elastic only) */
	/*************************************/

	Vector3r& us_elastic = phys->usElastic;
	us_elastic           = scg->rotate(us_elastic); // rotate vector
	Vector3r prevUs_el   = us_elastic;              // store previous elastic shear displacement (already rotated)
	us_elastic -= incidentVs * dt;                  // add shear increment

	/****************************************/
	/* SHEAR DISPLACEMENT (elastic+plastic) */
	/****************************************/

	Vector3r& us_total  = phys->usTotal;
	us_total            = scg->rotate(us_total); // rotate vector
	Vector3r prevUs_tot = us_total;              // store previous total shear displacement (already rotated)
	us_total -= incidentVs
	        * dt; // add shear increment NOTE: this vector is not passed into the failure criterion, hence it holds also the plastic part of the shear displacement

	bool noShearDamp = false; // bool to decide whether we need to account for shear damping dissipation or not

	/********************/
	/* MOHR-COULOMB law */
	/********************/
	phys->isSliding    = false;
	phys->shearViscous = Vector3r::Zero(); // reset so that during sliding, the previous values is not there
	Fn                 = phys->normalForce.norm();
	if (!includeAdhesion) {
		Real maxFs = Fn * phys->tangensOfFrictionAngle;
		if (shearElastic.squaredNorm() > maxFs * maxFs) {
			phys->isSliding = true;
			noShearDamp     = true; // no damping is added in the shear direction, hence no need to account for shear damping dissipation
			Real ratio      = maxFs / shearElastic.norm();
			shearElastic *= ratio;
			phys->shearForce = shearElastic; /*store only elastic shear displacement*/
			us_elastic *= ratio;
			if (calcEnergy) {
				frictionDissipation += (us_total - prevUs_tot).dot(shearElastic);
			}                                     // calculate energy dissipation due to sliding behavior
		} else if (useDamping) {                      // add current contact damping if we do not slide and if damping is requested
			phys->shearViscous = cs * incidentVs; // get shear viscous component
			phys->shearForce   = shearElastic - phys->shearViscous;
		} else if (!useDamping) {
			phys->shearForce = shearElastic;
		} // update the shear force at the elastic value if no damping is present and if we passed MC
	} else {  // Mohr-Coulomb formulation adpated due to the presence of adhesion (see Thornton, 1991).
		Real maxFs = phys->tangensOfFrictionAngle * (phys->adhesionForce + Fn); // adhesionForce already included in normalForce (above)
		if (shearElastic.squaredNorm() > maxFs * maxFs) {
			phys->isSliding = true;
			noShearDamp     = true; // no damping is added in the shear direction, hence no need to account for shear damping dissipation
			Real ratio      = maxFs / shearElastic.norm();
			shearElastic *= ratio;
			phys->shearForce = shearElastic; /*store only elastic shear displacement*/
			us_elastic *= ratio;
			if (calcEnergy) {
				frictionDissipation += (us_total - prevUs_tot).dot(shearElastic);
			}                                     // calculate energy dissipation due to sliding behavior
		} else if (useDamping) {                      // add current contact damping if we do not slide and if damping is requested
			phys->shearViscous = cs * incidentVs; // get shear viscous component
			phys->shearForce   = shearElastic - phys->shearViscous;
		} else if (!useDamping) {
			phys->shearForce = shearElastic;
		} // update the shear force at the elastic value if no damping is present and if we passed MC
	}

	/************************/
	/* SHEAR ELASTIC ENERGY */
	/************************/

	// NOTE: shear elastic energy calculation must come after the MC criterion, otherwise displacements and forces are not updated
	if (calcEnergy) {
		shearEnergy
		        += (us_elastic - prevUs_el)
		                   .dot((shearElastic + prev_FsElastic)
		                        / 2.); // NOTE: no additional energy if we perform sliding since us_elastic and prevUs_el will hold the same value (in fact us_elastic is only keeping the elastic part). We work out the area of the trapezium.
	}

	/**************************************************/
	/* VISCOUS DAMPING (energy term, shear direction) */
	/**************************************************/

	if (useDamping) { // get normal viscous component (the shear one is calculated inside Mohr-Coulomb criterion, see above)
		if (calcEnergy) {
			if (!noShearDamp) { shearDampDissip += phys->shearViscous.dot(incidentVs * dt); }
		} // calc energy dissipation due to viscous linear damping
	}

	/****************/
	/* APPLY FORCES */
	/****************/

	if (!scene->isPeriodic)
		applyForceAtContactPoint(-phys->normalForce - phys->shearForce, scg->contactPoint, id1, de1->se3.position, id2, de2->se3.position);
	else { // in scg we do not wrap particles positions, hence "applyForceAtContactPoint" cannot be used
		Vector3r force = -phys->normalForce - phys->shearForce;
		scene->forces.addForce(id1, force);
		scene->forces.addForce(id2, -force);
		scene->forces.addTorque(id1, (scg->radius1 - 0.5 * scg->penetrationDepth) * scg->normal.cross(force));
		scene->forces.addTorque(id2, (scg->radius2 - 0.5 * scg->penetrationDepth) * scg->normal.cross(force));
	}

	/********************************************/
	/* MOMENT CONTACT LAW */
	/********************************************/
	if (includeMoment) {
		// *** Bending ***//
		// new code to compute relative particle rotation (similar to the way the shear is computed)
		// use scg function to compute relAngVel
		Vector3r relAngVel = scg->getRelAngVel(de1, de2, dt);
		//Vector3r relAngVel = (b2->state->angVel-b1->state->angVel);
		Vector3r relAngVelBend = relAngVel - scg->normal.dot(relAngVel) * scg->normal; // keep only the bending part
		Vector3r relRot        = relAngVelBend * dt;                                   // relative rotation due to rolling behaviour
		// incremental formulation for the bending moment (as for the shear part)
		Vector3r& momentBend = phys->momentBend;
		momentBend           = scg->rotate(momentBend);        // rotate moment vector (updated)
		momentBend           = momentBend - phys->kr * relRot; // add incremental rolling to the rolling vector
		// ----------------------------------------------------------------------------------------
		// *** Torsion ***//
		Vector3r relAngVelTwist = scg->normal.dot(relAngVel) * scg->normal;
		Vector3r relRotTwist    = relAngVelTwist * dt; // component of relative rotation along n
		// incremental formulation for the torsional moment
		Vector3r& momentTwist = phys->momentTwist;
		momentTwist           = scg->rotate(momentTwist); // rotate moment vector (updated)
		momentTwist           = momentTwist - phys->ktw * relRotTwist;

#if 0
	// code to compute the relative particle rotation
	if (includeMoment){
		Real rMean = (scg->radius1+scg->radius2)/2.;
		// sliding motion
		Vector3r duS1 = scg->radius1*(phys->prevNormal-scg->normal);
		Vector3r duS2 = scg->radius2*(scg->normal-phys->prevNormal);
		// rolling motion
		Vector3r duR1 = scg->radius1*dt*b1->state->angVel.cross(scg->normal);
		Vector3r duR2 = -scg->radius2*dt*b2->state->angVel.cross(scg->normal);
		// relative position of the old contact point with respect to the new one
		Vector3r relPosC1 = duS1+duR1;
		Vector3r relPosC2 = duS2+duR2;

		Vector3r duR = (relPosC1+relPosC2)/2.; // incremental displacement vector (same radius is temporarily assumed)

		// check wheter rolling will be present, if not do nothing
		Vector3r x=scg->normal.cross(duR);
		Vector3r normdThetaR(Vector3r::Zero()); // initialize
		if(x.squaredNorm()==0) { /* no rolling */ }
		else {
				Vector3r normdThetaR = x/x.norm(); // moment unit vector
				phys->dThetaR = duR.norm()/rMean*normdThetaR;} // incremental rolling

		// incremental formulation for the bending moment (as for the shear part)
		Vector3r& momentBend = phys->momentBend;
		momentBend = scg->rotate(momentBend); // rotate moment vector
		momentBend = momentBend+phys->kr*phys->dThetaR; // add incremental rolling to the rolling vector FIXME: is the sign correct?
#endif

		// check plasticity condition (only bending part for the moment)
		Real MomentMax    = phys->maxBendPl * phys->normalForce.norm();
		Real scalarMoment = phys->momentBend.norm();
		if (MomentMax > 0) {
			if (scalarMoment > MomentMax) {
				Real ratio = MomentMax / scalarMoment; // to fix the moment to its yielding value
				phys->momentBend *= ratio;
			}
		}
		// apply moments
		Vector3r moment = phys->momentTwist + phys->momentBend;
		scene->forces.addTorque(id1, -moment);
		scene->forces.addTorque(id2, moment);
	}
	return true;
	// update variables
	//phys->prevNormal = scg->normal;
}

// The following code was moved from Ip2_FrictMat_FrictMat_MindlinCapillaryPhys.cpp

void Ip2_FrictMat_FrictMat_MindlinCapillaryPhys::go(
        const shared_ptr<Material>& b1 //FrictMat
        ,
        const shared_ptr<Material>& b2 // FrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return; // no updates of an already existing contact necessary

	shared_ptr<MindlinCapillaryPhys> contactPhysics(new MindlinCapillaryPhys());
	interaction->phys = contactPhysics;

	const auto mat1 = YADE_CAST<FrictMat*>(b1.get());
	const auto mat2 = YADE_CAST<FrictMat*>(b2.get());

	/* from interaction physics */
	const Real Ea = mat1->young;
	const Real Eb = mat2->young;
	const Real Va = mat1->poisson;
	const Real Vb = mat2->poisson;
	const Real fa = mat1->frictionAngle;
	const Real fb = mat2->frictionAngle;

	/* from interaction geometry */
	const auto scg = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	const Real Da  = scg->refR1 > 0 ? scg->refR1 : scg->refR2;
	const Real Db  = scg->refR2 > 0 ? scg->refR2 : scg->refR1;


	//Vector3r normal=scg->normal;  //The variable set but not used

	/* calculate stiffness coefficients */
	const Real Ga = Ea / (2 * (1 + Va));
	const Real Gb = Eb / (2 * (1 + Vb));
	const Real G  = 1.0 / ((2 - Va) / Ga + (2 - Vb) / Gb); //(Ga + Gb) / 2;                 // effective shear modulus
	//	const Real V             = (Va + Vb) / 2;                                                           // average of poisson's ratio
	const Real E             = Ea * Eb / ((1. - math::pow(Va, 2)) * Eb + (1. - math::pow(Vb, 2)) * Ea); // Young modulus
	const Real R             = Da * Db / (Da + Db);                                                     // equivalent radius
	const Real Rmean         = (Da + Db) / 2.;                                                          // mean radius
	const Real Kno           = 4. / 3. * E * sqrt(R);                                                   // coefficient for normal stiffness
	const Real Kso           = 8 * sqrt(R) * G;                                                         // coefficient for shear stiffness
	const Real frictionAngle = math::min(fa, fb);

	const Real Adhesion = 4. * Mathr::PI * R * gamma; // calculate adhesion force as predicted by DMT theory

	/* pass values calculated from above to MindlinCapillaryPhys */
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	//mindlinPhys->prevNormal = scg->normal; // used to compute relative rotation
	contactPhysics->kno           = Kno; // this is just a coeff
	contactPhysics->kso           = Kso; // this is just a coeff
	contactPhysics->adhesionForce = Adhesion;

	contactPhysics->kr        = krot;
	contactPhysics->ktw       = ktwist;
	contactPhysics->maxBendPl = eta * Rmean; // does this make sense? why do we take Rmean?

	/* compute viscous coefficients */
	if (en && betan) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinCapillaryPhys: only one of en, betan can be specified.");
	if (es && betas) throw std::invalid_argument("Ip2_FrictMat_FrictMat_MindlinCapillaryPhys: only one of es, betas can be specified.");

	// en or es specified
	if (en || es) {
		const Real h1  = -6.918798; // Fitting coefficients h_i from  Table 2 - Thornton et al. (2013).
		const Real h2  = -16.41105;
		const Real h3  = 146.8049;
		const Real h4  = -796.4559;
		const Real h5  = 2928.711;
		const Real h6  = -7206.864;
		const Real h7  = 11494.29;
		const Real h8  = -11342.18;
		const Real h9  = 6276.757;
		const Real h10 = -1489.915;

		// Consider same coefficient of restitution if only one is given (en or es)
		if (!en) { en = es; }
		if (!es) { es = en; }

		const Real En     = (*en)(mat1->id, mat2->id);
		const Real Es     = (*es)(mat1->id, mat2->id);
		const Real alphan = En
		        * (h1
		           + En * (h2 + En * (h3 + En * (h4 + En * (h5 + En * (h6 + En * (h7 + En * (h8 + En * (h9 + En * h10))))))))); // Eq. (B7) from Thornton et al. (2013)
		contactPhysics->betan = (En == 1.0) ? 0
		                                    : sqrt(1.0 / (1.0 - (math::pow(1.0 + En, 2)) * exp(alphan))
		                                           - 1.0); // Eq. (B6) from Thornton et al. (2013) - This is noted as 'gamma' in their paper

		// although Thornton (2015) considered betan=betas, here we use his formulae (B6) and (B7) allowing for betas to take a different value, based on the input es
		const Real alphas     = Es * (h1 + Es * (h2 + Es * (h3 + Es * (h4 + Es * (h5 + Es * (h6 + Es * (h7 + Es * (h8 + Es * (h9 + Es * h10)))))))));
		contactPhysics->betas = (Es == 1.0) ? 0 : sqrt(1.0 / (1.0 - (math::pow(1.0 + Es, 2)) * exp(alphas)) - 1.0);

		// betan/betas specified, use that value directly
	} else {
		contactPhysics->betan = betan ? (*betan)(mat1->id, mat2->id) : 0;
		contactPhysics->betas = betas ? (*betas)(mat1->id, mat2->id) : contactPhysics->betan;
	}
};

} // namespace yade
