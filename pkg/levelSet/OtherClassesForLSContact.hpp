/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/levelSet/LevelSet.hpp>

namespace yade {
class Bo1_LevelSet_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*) override;
	FUNCTOR1D(LevelSet);
	YADE_CLASS_BASE_DOC(Bo1_LevelSet_Aabb, BoundFunctor, "Creates/updates an :yref:`Aabb` of a :yref:`LevelSet`");
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Bo1_LevelSet_Aabb);


class MultiFrictPhys : public IPhys {
public:
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(MultiFrictPhys,IPhys,"A set of :yref:`FrictPhys` for describing the physical part of an interaction with multiple frictional contact points between two :yref:`LevelSet` bodies, as a set of :yref:`FrictPhys` items in :yref:`contacts<MultiFrictPhys.contacts>`. To combine with :yref:`MultiScGeom` and associated classes.",
	((vector< shared_ptr<FrictPhys> >,contacts,,,"The actual list of :yref:`FrictPhys` items corresponding to the different contact points."))
	((vector< int >,nodesIds,,,"The physics counterpart of :yref:`MultiScGeom.nodesIds` (both should be equal by design).")) // do we need both ?
	((Real,kn,0,,"Mother value of :yref:`FrictPhys.kn` that will apply to each contact point."))
	((Real,ks,0,,"Mother value of :yref:`FrictPhys.ks` that will apply to each contact point."))
	((Real,frictAngle,0,,"Mother value of atan(:yref:`FrictPhys.tangensOfFrictionAngle`) in radians that will apply to each contact point."))
	,
	createIndex(); // this class will enter InteractionLoop dispatch, we need a create_index() here, and a REGISTER_*_INDEX below (https://yade-dem.org/doc/prog.html#indexing-dispatch-types)
	);
	// clang-format on
	REGISTER_CLASS_INDEX(MultiFrictPhys, IPhys); // see createIndex() remark
};
REGISTER_SERIALIZABLE(MultiFrictPhys);

class Ip2_FrictMat_FrictMat_MultiFrictPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	FUNCTOR2D(FrictMat, FrictMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ip2_FrictMat_FrictMat_MultiFrictPhys,IPhysFunctor,"Create a :yref:`MultiFrictPhys` from two :yref:`FrictMats<FrictMat>`. Mother contact stiffnesses (:yref:`MultiFrictPhys.kn` and :yref:`MultiFrictPhys.ks`) are directly assigned from below attributes, independent of FrictMat properties. Global friction angle (:yref:`MultiFrictPhys.frictAngle`) is taken as the minimum of the 2 material friction angles (:yref:`FrictMat.frictionAngle`).",
		((Real,kn,0,,"Chosen value for :yref:`MultiFrictPhys.kn`"))
		((Real,ks,0,,"Chosen value for :yref:`MultiFrictPhys.ks`"))
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_MultiFrictPhys);

} // namespace yade
#endif // YADE_LS_DEM
