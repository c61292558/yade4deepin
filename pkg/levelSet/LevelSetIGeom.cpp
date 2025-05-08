/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <pkg/levelSet/LevelSetIGeom.hpp>

namespace yade {
YADE_PLUGIN((MultiScGeom));

bool MultiScGeom::hasNode(int nodeIdx) { return (iteratorToNode(nodeIdx) != nodesIds.end()); }

std::vector<int>::iterator MultiScGeom::iteratorToNode(int nodeIdx)
{
	return (std::find(nodesIds.begin(), nodesIds.end(), nodeIdx) // no problem even if nodesIds is empty)
	);
}

} // namespace yade
#endif //YADE_LS_DEM
