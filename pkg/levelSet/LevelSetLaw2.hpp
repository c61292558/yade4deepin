/*************************************************************************
*  2021-24 jerome.duriez@inrae.fr                                        *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <pkg/dem/ElasticContactLaw.hpp>
#include <pkg/levelSet/LevelSetIGeom.hpp>            // for MultiScGeom
#include <pkg/levelSet/OtherClassesForLSContact.hpp> // for MultiFrictPhys

namespace yade { // Cannot have #include directive inside.

class Law2_MultiScGeom_MultiFrictPhys_CundallStrack : public Law2_ScGeom_FrictPhys_CundallStrack {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
  YADE_CLASS_BASE_DOC(Law2_MultiScGeom_MultiFrictPhys_CundallStrack,Law2_ScGeom_FrictPhys_CundallStrack,"Applies :yref:`Law2_ScGeom_FrictPhys_CundallStrack` at each contact point of a (yref:`MultiScGeom`;yref:`MultiFrictPhys`) contact [Duriez2023]_."//,TODO: check (and correct ?) behavior of Python .def-ed functions being inherited from mother class
	);
	// clang-format on
	FUNCTOR2D(MultiScGeom, MultiFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_MultiScGeom_MultiFrictPhys_CundallStrack);

} // namespace yade
#endif //YADE_LS_DEM
