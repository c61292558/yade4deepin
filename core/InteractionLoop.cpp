#include "InteractionLoop.hpp"

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.


YADE_PLUGIN((InteractionLoop));
CREATE_LOGGER(InteractionLoop);

void InteractionLoop::pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& /*d*/)
{
	if (boost::python::len(t) == 0) return; // nothing to do
	if (boost::python::len(t) != 3) throw invalid_argument("Exactly 3 lists of functors must be given");
	// parse custom arguments (3 lists) and do in-place modification of args
	using vecGeom = std::vector<shared_ptr<IGeomFunctor>>;
	using vecPhys = std::vector<shared_ptr<IPhysFunctor>>;
	using vecLaw  = std::vector<shared_ptr<LawFunctor>>;
	vecGeom vg    = boost::python::extract<vecGeom>(t[0])();
	vecPhys vp    = boost::python::extract<vecPhys>(t[1])();
	vecLaw  vl    = boost::python::extract<vecLaw>(t[2])();
	for (const auto& gf : vg)
		this->geomDispatcher->add(gf);
	for (const auto& pf : vp)
		this->physDispatcher->add(pf);
	for (const auto& cf : vl)
		this->lawDispatcher->add(cf);
	t = boost::python::tuple(); // empty the args; not sure if this is OK, as there is some refcounting in raw_constructor code
}

void InteractionLoop::action()
{
	// update Scene* of the dispatchers
	lawDispatcher->scene  = scene;
	physDispatcher->scene = scene;
	geomDispatcher->scene = scene;

	// ask dispatchers to update Scene* of their functors
	geomDispatcher->updateScenePtr();
	physDispatcher->updateScenePtr();
	lawDispatcher->updateScenePtr();

#ifdef YADE_MPI
	const Body::id_t& subDIdx = scene->subdomain;
#endif

	/*
		initialize callbacks; they return pointer (used only in this timestep) to the function to be called
		returning NULL deactivates the callback in this timestep
	*/
	// pair of callback object and pointer to the function to be called
	vector<IntrCallback::FuncPtr> callbackPtrs;
	for (const auto& cb : callbacks) {
		cb->scene = scene;
		callbackPtrs.push_back(cb->stepInit());
	}
	assert(callbackPtrs.size() == callbacks.size());
	const size_t callbacksSize = callbacks.size();

	// cache transformed cell size
	Matrix3r cellHsize = Matrix3r::Zero();
	if (scene->isPeriodic) { cellHsize = scene->cell->hSize; }

	// force removal of interactions that were not encountered by the collider
	// (only for some kinds of colliders; see comment for InteractionContainer::iterColliderLastRun)
	const bool removeUnseenIntrs = (scene->interactions->iterColliderLastRun >= 0 && scene->interactions->iterColliderLastRun == scene->iter);

	const long size = scene->interactions->size();

	vector<shared_ptr<Interaction>>* interactions; //a pointer to an interaction vector.
	if (loopOnSortedInteractions) {
		scene->interactions->updateSortedIntrs();           //sort sortedIntrs, this is VERY SLOW !
		interactions = &(scene->interactions->sortedIntrs); //set the pointer to the address of the sorted version of the vector
	} else
		interactions = &(scene->interactions
		                         ->linIntrs); //set the pointer to the address of the unsorted version of the vector (original version, normal behavior)
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& I = (*interactions)[i];
		if (removeUnseenIntrs && !I->isReal() && I->iterLastSeen < scene->iter) {
			eraseAfterLoop(I->getId1(), I->getId2());
			continue;
		}
		const shared_ptr<Body>& b1_ = Body::byId(I->getId1(), scene);
		const shared_ptr<Body>& b2_ = Body::byId(I->getId2(), scene);

		if (!b1_ || !b2_) {
			// 			LOG_DEBUG("Body #"<<(b1_?I->getId2():I->getId1())<<" vanished, erasing intr #"<<I->getId1()<<"+#"<<I->getId2()<<"!");
			scene->interactions->requestErase(I);
			continue;
		}

#ifdef YADE_MPI
		// Skip interactions between remote bodies, and reset them so we don't keep deprecated data
		if (subDIdx != b1_->subdomain and subDIdx != b2_->subdomain) {
			scene->interactions->requestErase(I->getId1(), I->getId2());
			continue;
		}
#endif

		// Skip interaction with clumps
		if (b1_->isClump() || b2_->isClump()) { continue; }

		// we know there is no geometry functor already, take the short path
		if (!I->functorCache.geomExists) {
			assert(!I->isReal());
			continue;
		}

		// no interaction geometry for either of bodies; no interaction possible
		if (!b1_->shape || !b2_->shape) {
			assert(!I->isReal());
			continue;
		}

		bool swap = false;
		// IGeomDispatcher
		if (!I->functorCache.geom) {
			I->functorCache.geom = geomDispatcher->getFunctor2D(b1_->shape, b2_->shape, swap);
			// returns NULL ptr if no functor exists; remember that and shortcut
			if (!I->functorCache.geom) {
				I->functorCache.geomExists = false;
				continue;
			}
		}
		// arguments for the geom functor are in the reverse order (dispatcher would normally call goReverse).
		// we don't remember the fact that is reverse, so we swap bodies within the interaction
		// and can call go in all cases
		if (swap) { I->swapOrder(); }
		// body pointers must be updated, in case we swapped
		const shared_ptr<Body>& b1 = swap ? b2_ : b1_;
		const shared_ptr<Body>& b2 = swap ? b1_ : b2_;

		assert(I->functorCache.geom);

		bool wasReal = I->isReal();
		bool geomCreated;
		if (!scene->isPeriodic) {
			geomCreated = I->functorCache.geom->go(b1->shape, b2->shape, *b1->state, *b2->state, Vector3r::Zero(), /*force*/ false, I);
		} else {
			// handle periodicity
			Vector3r shift2 = cellHsize * I->cellDist.cast<Real>();
			// in sheared cell, apply shear on the mutual position as well
			geomCreated = I->functorCache.geom->go(b1->shape, b2->shape, *b1->state, *b2->state, shift2, /*force*/ false, I);
		}
		if (!geomCreated) {
			if (wasReal) LOG_WARN("IGeomFunctor returned false on existing interaction!");
			if (wasReal)
				scene->interactions->requestErase(I); // fully created interaction without geometry is reset and perhaps erased in the next step
			continue;                                     // in any case don't care about this one anymore
		}

		// IPhysDispatcher
		if (!I->functorCache.phys) {
			I->functorCache.phys = physDispatcher->getFunctor2D(b1->material, b2->material, swap);
			assert(!swap); // InteractionPhysicsEngineUnits are symmetric
		}

		if (!I->functorCache.phys) {
			throw std::runtime_error(
			        "Undefined or ambiguous IPhys dispatch for types " + b1->material->getClassName() + " and " + b2->material->getClassName()
			        + ".");
		}
		I->functorCache.phys->go(b1->material, b2->material, I);
		assert(I->phys);

		if (!wasReal) I->iterMadeReal = scene->iter; // mark the interaction as created right now

		// LawDispatcher
		// populating constLaw cache must be done after geom and physics dispatchers have been called, since otherwise the interaction
		// would not have geom and phys yet.
		if (!I->functorCache.constLaw) {
			I->functorCache.constLaw = lawDispatcher->getFunctor2D(I->geom, I->phys, swap);
			if (!I->functorCache.constLaw) {
				LOG_FATAL(
				        "None of given Law2 functors can handle interaction #"
				        << I->getId1() << "+" << I->getId2() << ", types geom:" << I->geom->getClassName() << "=" << I->geom->getClassIndex()
				        << " and phys:" << I->phys->getClassName() << "=" << I->phys->getClassIndex()
				        << " (LawDispatcher::getFunctor2D returned empty functor)");
				exit(1);
			}
			assert(!swap); // reverse call would make no sense, as the arguments are of different types
		}
		assert(I->functorCache.constLaw);

		//If the functor return false, the interaction is reset
		if (!I->functorCache.constLaw->go(I->geom, I->phys, I.get())) scene->interactions->requestErase(I);

		// process callbacks for this interaction
		// 		Note: the following condition is algorithmicaly safe, however a possible use of callbacks is to do something special when interactions are deleted, which is impossible if we skip them. The test should be commented out
		if (!I->isReal()) continue; // it is possible that Law2_ functor called requestErase, hence this check
		for (size_t j = 0; j < callbacksSize; j++) {
			if (callbackPtrs[j] != NULL) (*(callbackPtrs[j]))(callbacks[j].get(), I.get());
		}
	}
}

shared_ptr<Interaction> InteractionLoop::createExplicitInteraction(Body::id_t id1, Body::id_t id2, bool force, bool virtualI)
{
	IGeomDispatcher*        geomMeta = NULL;
	IPhysDispatcher*        physMeta = NULL;
	shared_ptr<Scene>       rb       = Omega::instance().getScene();
	shared_ptr<Interaction> i        = rb->interactions->find(Body::id_t(id1), Body::id_t(id2));
	if (i) {
		if (i->isReal())
			throw runtime_error(
			        string("Interaction #") + boost::lexical_cast<string>(id1) + "+#" + boost::lexical_cast<string>(id2) + " already exists.");
		else
			rb->interactions->erase(id1, id2, i->linIx);
	}
	shared_ptr<Body> b1 = Body::byId(id1, rb), b2 = Body::byId(id2, rb);
	if (!b1) throw runtime_error(("No body #" + boost::lexical_cast<string>(id1)).c_str());
	if (!b2) throw runtime_error(("No body #" + boost::lexical_cast<string>(id2)).c_str());

	if (not virtualI) { // normal case, create a statefull interaction or nothing
		FOREACH(const shared_ptr<Engine>& e, rb->engines)
		{
			if (!geomMeta) {
				geomMeta = dynamic_cast<IGeomDispatcher*>(e.get());
				if (geomMeta) continue;
			}
			if (!physMeta) {
				physMeta = dynamic_cast<IPhysDispatcher*>(e.get());
				if (physMeta) continue;
			}
			InteractionLoop* id(dynamic_cast<InteractionLoop*>(e.get()));
			if (id) {
				geomMeta = id->geomDispatcher.get();
				physMeta = id->physDispatcher.get();
			}
			if (geomMeta && physMeta) { break; }
		}
		if (!geomMeta) throw runtime_error("No IGeomDispatcher in engines or inside InteractionLoop.");
		if (!physMeta) throw runtime_error("No IPhysDispatcher in engines or inside InteractionLoop.");
		i = geomMeta->explicitAction(b1, b2, /*force*/ force);
		assert(force && i);
		if (!i) return i;
		physMeta->explicitAction(b1->material, b2->material, i);
		i->iterMadeReal = rb->iter;
	} else
		i = shared_ptr<Interaction>(new Interaction(id1, id2));
	rb->interactions->insert(i);
	return i;
}

} // namespace yade
