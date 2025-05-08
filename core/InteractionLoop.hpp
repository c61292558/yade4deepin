// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include <core/Callbacks.hpp>
#include <core/Dispatching.hpp>
#include <core/GlobalEngine.hpp>

#ifdef USE_TIMING_DELTAS
#define TIMING_DELTAS_CHECKPOINT(cpt) timingDeltas->checkpoint(cpt)
#define TIMING_DELTAS_START() timingDeltas->start()
#else
#define TIMING_DELTAS_CHECKPOINT(cpt)
#define TIMING_DELTAS_START()
#endif

namespace yade { // Cannot have #include directive inside.


class InteractionLoop : public GlobalEngine {
	bool alreadyWarnedNoCollider;
	using idPair = std::pair<Body::id_t, Body::id_t>;
// store interactions that should be deleted after loop in action, not later
#ifdef YADE_OPENMP
	std::vector<std::list<idPair>> eraseAfterLoopIds;
	void                           eraseAfterLoop(Body::id_t id1, Body::id_t id2) { eraseAfterLoopIds[omp_get_thread_num()].push_back(idPair(id1, id2)); }
#else
	std::list<idPair> eraseAfterLoopIds;
	void              eraseAfterLoop(Body::id_t id1, Body::id_t id2) { eraseAfterLoopIds.push_back(idPair(id1, id2)); }
	//! create transientInteraction between 2 bodies, using existing Dispatcher in Omega
#endif
public:
	void                           pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& d) override;
	static shared_ptr<Interaction> createExplicitInteraction(Body::id_t id1, Body::id_t id2, bool force, bool virtualI);
	void                           action() override;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(InteractionLoop,GlobalEngine,"Unified dispatcher for handling interaction loop at every step, for parallel performance reasons.\n\n.. admonition:: Special constructor\n\n\tConstructs from 3 lists of :yref:`Ig2<IGeomFunctor>`, :yref:`Ip2<IPhysFunctor>`, :yref:`Law2<LawFunctor>` functors respectively; they will be passed to internal dispatchers, which you might retrieve as :yref:`geomDispatcher<InteractionLoop.geomDispatcher>`, :yref:`physDispatcher<InteractionLoop.physDispatcher>`, :yref:`lawDispatcher<InteractionLoop.lawDispatcher>` respectively.",
			((shared_ptr<IGeomDispatcher>,geomDispatcher,new IGeomDispatcher,Attr::readonly,":yref:`IGeomDispatcher` object that is used for dispatch."))
			((shared_ptr<IPhysDispatcher>,physDispatcher,new IPhysDispatcher,Attr::readonly,":yref:`IPhysDispatcher` object used for dispatch."))
			((shared_ptr<LawDispatcher>,lawDispatcher,new LawDispatcher,Attr::readonly,":yref:`LawDispatcher` object used for dispatch."))
			((vector<shared_ptr<IntrCallback> >,callbacks,,,":yref:`Callbacks<IntrCallback>` which will be called for every :yref:`Interaction`, if activated."))
			((bool, loopOnSortedInteractions, false,,"If true, the main interaction loop will occur on a sorted list of interactions. This is SLOW but useful to workaround floating point force addition non reproducibility when debugging parallel implementations of yade."))
			,
			/*ctor*/ alreadyWarnedNoCollider=false;
				#ifdef YADE_OPENMP
					eraseAfterLoopIds.resize(omp_get_max_threads());
				#endif
			,
			/*py*/
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(InteractionLoop);

} // namespace yade
