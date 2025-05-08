/*************************************************************************
*  2021-2024 jerome.duriez@inrae.fr                                      *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/IGeom.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade {

class MultiScGeom : public IGeom {
public:
	// hasNode and iteratorToNode below could certainly be const-declared but this may require a better handling of int / const int argument, maybe with respect to find in iteratorToNode which would expect a const int as last argument
	bool                       hasNode(int);        // to check whether some node (idx) is in nodesIds. Could there be sthg as simple as Python "in" ?
	std::vector<int>::iterator iteratorToNode(int); // for locating some node (idx) in nodesIds
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(MultiScGeom,IGeom,"A set of :yref:`ScGeom` for describing the kinematics of an interaction with multiple contact points between two :yref:`LevelSet` bodies, as a set of :yref:`ScGeom` items in :yref:`contacts<MultiScGeom.contacts>`. To combine with :yref:`MultiFrictPhys` and associated classes.",
	((vector< shared_ptr<ScGeom> >,contacts,,,"The actual list of :yref:`ScGeom` items corresponding to the different contact points."))
	((vector< int >,nodesIds,,,"List of :yref:`surface nodes<LevelSet.surfaceNodes>` (on id1 if that body is smaller -- or equal -- in volume, or id2 otherwise) making :yref:`contacts<MultiScGeom.contacts>`. Contact point for a node of index nodesIds[i] has kinematic properties stored in contacts[i]. Should be equal to :yref:`MultiFrictPhys.nodesIds` by design")) // do we need both ?"
	, // existing Law2::go() will have to apply on a reference but contacts itself can not be a vector of references: https://stackoverflow.com/questions/922360/why-cant-i-make-a-vector-of-references
	createIndex(); // this class will enter InteractionLoop dispatch, we need a create_index() here, and a REGISTER_*_INDEX below (https://yade-dem.org/doc/prog.html#indexing-dispatch-types)
	);
	// clang-format on
	DECLARE_LOGGER;
	REGISTER_CLASS_INDEX(MultiScGeom, IGeom); // see createIndex() remark
};
REGISTER_SERIALIZABLE(MultiScGeom);

} // namespace yade
#endif // YADE_LS_DEM
