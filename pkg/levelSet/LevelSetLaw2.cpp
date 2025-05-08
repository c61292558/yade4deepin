/*************************************************************************
*  2021-24 jerome.duriez@inrae.fr                                        *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <pkg/levelSet/LevelSetLaw2.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((Law2_MultiScGeom_MultiFrictPhys_CundallStrack));
CREATE_LOGGER(Law2_MultiScGeom_MultiFrictPhys_CundallStrack);

bool Law2_MultiScGeom_MultiFrictPhys_CundallStrack::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	shared_ptr<MultiScGeom>    igMulti = YADE_PTR_CAST<MultiScGeom>(ig);
	shared_ptr<MultiFrictPhys> ipMulti = YADE_PTR_CAST<MultiFrictPhys>(ip);
	shared_ptr<IGeom>          igOne(
                new IGeom); // keep it as IGeom ! Otherwise (if directly defined as ScGeom) a temporary will have to be created when calling CS::go below, and it will not work with the expected reference type
	shared_ptr<IPhys> ipOne(new IPhys);
	bool              retVal(false);
	LOG_DEBUG("Will loop over " << igMulti->contacts.size() << " contacts looking at igMulti, vs " << ipMulti->contacts.size() << " in ipMulti");
	for (unsigned int idx = 0; idx < igMulti->contacts.size(); idx++) {
		LOG_TRACE("Looping over contact " << idx << " out of " << igMulti->contacts.size());
		igOne  = igMulti->contacts[idx];
		ipOne  = ipMulti->contacts[idx];
		retVal = Law2_ScGeom_FrictPhys_CundallStrack::go(igOne, ipOne, contact) || retVal;
		// NB: in the above, value of Law2_MultiScGeom_MultiFrictPhys_CundallStrack.sphericalBodies is indeed taken into account
	}
	return retVal;
}

} // namespace yade
#endif //YADE_LS_DEM
