#include <core/Dispatching.hpp>
#include <core/InteractionLoop.hpp>
#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.


YADE_PLUGIN((BoundFunctor)(IGeomFunctor)(IPhysFunctor)(LawFunctor)(BoundDispatcher)(IGeomDispatcher)(IPhysDispatcher)(LawDispatcher));
BoundFunctor::~BoundFunctor() {};
IGeomFunctor::~IGeomFunctor() {};
IPhysFunctor::~IPhysFunctor() {};
LawFunctor::~LawFunctor() {};


/********************************************************************
                      BoundDispatcher
*********************************************************************/

CREATE_LOGGER(BoundDispatcher);
void BoundDispatcher::action()
{
	updateScenePtr();
	shared_ptr<BodyContainer>& bodies   = scene->bodies;
	const bool                 redirect = bodies->useRedirection;
	if (redirect) bodies->updateRealBodies();
	const long numBodies = redirect ? (long)bodies->realBodies.size() : (long)bodies->size();
#ifdef YADE_MPI
	Body::id_t subdomainId = 0;
#endif
#ifdef YADE_OPENMP
#pragma omp parallel for num_threads(ompThreads > 0 ? min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
	for (int id = 0; id < numBodies; id++) {
		if (not redirect and not bodies->exists(id)) continue; // don't delete this check  - Janek
		const shared_ptr<Body>& b = (*bodies)[redirect ? bodies->realBodies[id] : id];
		processBody(b);
#ifndef YADE_MPI
	}
#else //when all ordinary bodies have been processed, we will evaluate subdomain's min and max
		if (b->getIsSubdomain() and b->subdomain == scene->subdomain)
			subdomainId = b->getId(); //subdomain bounds need all other bodies to have updated bounds, hence we keep it for after this loop
	}
	if (subdomainId != 0) {
		YADE_PTR_CAST<Subdomain>((*bodies)[subdomainId]->shape)->setMinMax();
		processBody((*bodies)[subdomainId]);
	}
#endif
}

void BoundDispatcher::processBody(const shared_ptr<Body>& b)
{
	shared_ptr<Shape>& shape = b->shape;
	if (!b->isBounded() || !shape) return;
#ifdef BV_FUNCTOR_CACHE
	if (!shape->boundFunctor) {
		shape->boundFunctor = this->getFunctor1D(shape);
		if (!shape->boundFunctor) return;
	}
	shape->boundFunctor->go(shape, b->bound, b->state->se3, b.get());
#else
	operator()(shape, b->bound, b->state->se3, b.get());
#endif
	if (!b->bound) return; // the functor did not create new bound
	Real& sweepLength = b->bound->sweepLength;
	if (targetInterv > 0 and scene->iter > b->bound->lastUpdateIter) { //at iteration zero checking displacement makes no sense
		Vector3r disp = b->state->pos - b->bound->refPos;
		Real     dist = max(math::abs(disp[0]), max(math::abs(disp[1]), math::abs(disp[2])));
		if (dist or b->state->vel != Vector3r::Zero()) {
			Real newLength = dist * targetInterv / (scene->iter - b->bound->lastUpdateIter);
			newLength      = max(0.9 * sweepLength, newLength); //don't decrease size too fast to prevent time consuming oscillations
			sweepLength    = max(minSweepDistFactor * sweepDist, min(newLength, sweepDist));
		} else {
			sweepLength = 0;
		}
	} else
		sweepLength = sweepDist;
#ifdef YADE_MPI
	if (b->getIsSubdomain()) sweepLength = 0;
	// skip fluid mesh bounding box from being extended
	if (b->getIsFluidDomainBbox()) sweepLength = 0;
#endif
	b->bound->refPos         = b->state->pos;
	b->bound->lastUpdateIter = scene->iter;
	if (sweepLength > 0) {
		Aabb* aabb = YADE_CAST<Aabb*>(b->bound.get());
		aabb->min -= Vector3r(sweepLength, sweepLength, sweepLength);
		aabb->max += Vector3r(sweepLength, sweepLength, sweepLength);
	}
}


/********************************************************************
                      IGeomDispatcher
*********************************************************************/

CREATE_LOGGER(IGeomDispatcher);

shared_ptr<Interaction> IGeomDispatcher::explicitAction(const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool force)
{
	scene             = Omega::instance().getScene().get(); // to make sure if called from outside of the loop
	Vector3i cellDist = Vector3i::Zero();
	if (scene->isPeriodic) {
		for (int i = 0; i < 3; i++)
			cellDist[i] = -(int)((b2->state->pos[i] - b1->state->pos[i]) / scene->cell->getSize()[i] + .5);
	}
	Vector3r shift2 = scene->cell->hSize * cellDist.cast<Real>();
	updateScenePtr();

	assert(b1->shape && b2->shape);
	shared_ptr<Interaction> I(new Interaction(b1->getId(), b2->getId()));
	I->cellDist = cellDist;
	// FIXME: this code is more or less duplicated from InteractionLoop :-(
	bool swap            = false;
	I->functorCache.geom = getFunctor2D(b1->shape, b2->shape, swap);
	if (!I->functorCache.geom)
		throw invalid_argument(
		        "IGeomDispatcher::explicitAction could not dispatch for given types (" + b1->shape->getClassName() + "," + b2->shape->getClassName()
		        + ").");
	if (swap) { I->swapOrder(); }
	const shared_ptr<Body>& b1Swp = Body::byId(I->getId1(), scene);
	const shared_ptr<Body>& b2Swp = Body::byId(I->getId2(), scene);
	bool                    succ  = I->functorCache.geom->go(b1Swp->shape, b2Swp->shape, *b1Swp->state, *b2Swp->state, shift2, /*force*/ true, I);
	if (!succ and force)
		throw logic_error(
		        "Functor " + I->functorCache.geom->getClassName() + "::go returned false, even if asked to force IGeom creation. Please report bug.");
	return I;
}

void IGeomDispatcher::action()
{
	updateScenePtr();

	shared_ptr<BodyContainer>& bodies = scene->bodies;
	const bool                 isPeriodic(scene->isPeriodic);
	Matrix3r                   cellHsize;
	if (isPeriodic) cellHsize = scene->cell->hSize;
	bool removeUnseenIntrs = (scene->interactions->iterColliderLastRun >= 0 && scene->interactions->iterColliderLastRun == scene->iter);
#ifdef YADE_OPENMP
	const long size = scene->interactions->size();
#pragma omp parallel for
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& I = (*scene->interactions)[i];
#else
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
#endif
		if (removeUnseenIntrs && !I->isReal() && I->iterLastSeen < scene->iter) {
			scene->interactions->requestErase(I);
			continue;
		}
		const shared_ptr<Body>& b1 = (*bodies)[I->getId1()];
		const shared_ptr<Body>& b2 = (*bodies)[I->getId2()];

		if (!b1 || !b2) {
			// This code is duplicated in Dispatching.cpp:123 and InteractionLoop.cpp:73
			LOG_DEBUG("Body #" << (b1 ? I->getId2() : I->getId1()) << " vanished, erasing intr #" << I->getId1() << "+#" << I->getId2() << "!");
			scene->interactions->requestErase(I);
			continue;
		}

		bool wasReal = I->isReal();
		if (!b1->shape || !b2->shape) {
			assert(!wasReal);
			continue;
		} // some bodies do not have shape
		bool geomCreated;
		//const bool forceFalse=false;
		if (!isPeriodic) {
			geomCreated = I->functorCache.geom->go(b1->shape, b2->shape, *b1->state, *b2->state, Vector3r::Zero(), /*force*/ false, I);
		} else {
			Vector3r shift2 = cellHsize * I->cellDist.cast<Real>();
			geomCreated     = I->functorCache.geom->go(b1->shape, b2->shape, *b1->state, *b2->state, shift2, /*force*/ false, I);
		}
		// reset && erase interaction that existed but now has no geometry anymore
		if (wasReal && !geomCreated) { scene->interactions->requestErase(I); }
	}
}

/********************************************************************
                      IPhysDispatcher
*********************************************************************/


void IPhysDispatcher::explicitAction(shared_ptr<Material>& pp1, shared_ptr<Material>& pp2, shared_ptr<Interaction>& I)
{
	updateScenePtr();
	if (!I->geom) throw invalid_argument(string(__FILE__) + ": explicitAction received interaction without geom.");
	if (!I->functorCache.phys) {
		bool dummy;
		I->functorCache.phys = getFunctor2D(pp1, pp2, dummy);
		if (!I->functorCache.phys)
			throw invalid_argument(
			        "IPhysDispatcher::explicitAction did not find a suitable dispatch for types " + pp1->getClassName() + " and "
			        + pp2->getClassName());
		I->functorCache.phys->go(pp1, pp2, I);
	}
}

void IPhysDispatcher::action()
{
	updateScenePtr();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
#ifdef YADE_OPENMP
	const long size = scene->interactions->size();
#pragma omp parallel for
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& interaction = (*scene->interactions)[i];
#else
	FOREACH(const shared_ptr<Interaction>& interaction, *scene->interactions)
	{
#endif
		if (interaction->geom) {
			shared_ptr<Body>& b1      = (*bodies)[interaction->getId1()];
			shared_ptr<Body>& b2      = (*bodies)[interaction->getId2()];
			bool              hadPhys = (interaction->phys.get() != 0);
			                  operator()(b1->material, b2->material, interaction);
			assert(interaction->phys);
			if (!hadPhys) interaction->iterMadeReal = scene->iter;
		}
	}
}


/********************************************************************
                      LawDispatcher
*********************************************************************/

CREATE_LOGGER(LawDispatcher);
void LawDispatcher::action()
{
	updateScenePtr();
#ifdef YADE_OPENMP
	const long size = scene->interactions->size();
#pragma omp parallel for
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& I = (*scene->interactions)[i];
#else
	FOREACH(shared_ptr<Interaction> I, *scene->interactions)
	{
#endif
		if (I->isReal()) {
			assert(I->geom);
			assert(I->phys);
			operator()(I->geom, I->phys, I.get());
			if (!I->isReal() && I->isFresh(scene))
				LOG_ERROR("Law functor deleted interaction that was just created. Please report bug: either this message is spurious, or the "
				          "functor (or something else) is buggy.");
		}
	}
}

} // namespace yade
