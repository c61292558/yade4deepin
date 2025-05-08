/****************************************************************************
*  2023 DLH van der Haven, dannyvdhaven@gmail.com, University of Cambridge  *
*                                                                           *
*  For details, please see van der Haven et al., "A physically consistent   *
*   Discrete Element Method for arbitrary shapes using Volume-interacting   *
*   Level Sets", Comput. Methods Appl. Mech. Engrg., 414 (116165):1-21      *
*   https://doi.org/10.1016/j.cma.2023.116165                               *
*  This project has been financed by Novo Nordisk A/S (Bagsværd, Denmark).  *
*                                                                           *
*  This program is licensed under GNU GPLv2, see file LICENSE for details.  *
****************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/levelSet/LevelSet.hpp>
#include <pkg/levelSet/VolumeGeom.hpp>

namespace yade {

class Ig2_LevelSet_LevelSet_VolumeGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&,
	        const shared_ptr<Shape>&,
	        const State&,
	        const State&,
	        const Vector3r&,
	        const bool&,
	        const shared_ptr<Interaction>&)
	        override; // Reminder: method signature is imposed by InteractionLoop.cpp and also somewhat inherited from template class FunctorWrapper.
	bool
	goReverse(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override
	{
		LOG_ERROR("We ended up calling goReverse. How is this possible for symmetric IgFunctor ? Anyway, we now have to code something.");
		return false;
	};
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ig2_LevelSet_LevelSet_VolumeGeom,IGeomFunctor,R"""(Creates or updates a :yref:`VolumeGeom` instance representing the contact of two :yref:`LevelSet` bodies of arbitrary shape.	An algorithm is used that recursively evaluates the signed distance function $\phi$ (:yref:`LevelSet.distField`) at increasingly finer mesh sizes to compute the overlap volume $V$. Surface nodes are obsolete if this functor is used. Denoting $u_n$ as the :yref:`overlap<VolumeGeom.penetrationVolume>`, $\vec{C}$ the :yref:`contact point<VolumeGeom.contactPoint>` and $\vec{n}$ the contact :yref:`normal<VolumeGeom.normal>`, we have:

* $u_n = V_n = \sum_i V_i(\vec{x}_i)$
* $\vec{n} = \frac{ \sum_i V_i \vec{\nabla} \phi_1(\vec{x_i}) - \sum_i V_i \vec{\nabla} \phi_2(\vec{x_i}) }{2V_n}$  chosen to be oriented from :yref:`1<Interaction.id1>` to :yref:`2<Interaction.id2>`
* $\vec{C} = \frac{1}{V_n} \sum_i \vec{x}_i * V_i(\vec{x}_i)$

.. note:: Because this functor expresses the particle overlap $u_n$ as a volume, care needs to be taken that only volume-based contact laws are used. Current contact laws in YADE mainly use the overlap distance to determine the force.
)""",
	((uint,nRefineOctree,5,,"The number of refinements performed by the Octree algorithm used to compute the overlap volume between two particles. Default is 5. Note: (nr of layers, effective nr of integration elements): (1,1), (2,8), (3,64), (4,512), (5,4096), (n,8^(n-1))."))
	((Real,smearCoeffOctree,1.0,,"Smearing coefficient for the smeared Heaviside step function in the overlap volume integration. The transition width, or smearing width, is equal to half the diagonal of the smallest integration cell divided by the smearing coefficient."))
	((bool,useAABE,false,,"If true, use the provided (locally) axis-aligned bounding ellipsoid (AABE) to reduce the potential overlap volume between the particles. Increases accuracy of the Octree algrithm because the smallest integration cells will be smaller."))
	,,
	);
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(LevelSet, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(LevelSet, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_LevelSet_LevelSet_VolumeGeom);


class Ig2_Wall_LevelSet_VolumeGeom : public IGeomFunctor {
public:
	bool go(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override;
	bool goReverse(
	        const shared_ptr<Shape>&       cm1,
	        const shared_ptr<Shape>&       cm2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& c) override
	{
		c->swapOrder();
		return go(cm2, cm1, state2, state1, -shift2, force, c);
	};
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ig2_Wall_LevelSet_VolumeGeom,IGeomFunctor,"Creates or updates a :yref:`VolumeGeom` instance representing the intersection of one :yref:`LevelSet` body with one :yref:`Wall` body, where overlap is chosen to occur on the opposite wall side than the LevelSet body's center. :yref:`Contact normal<VolumeGeom.normal>` is given by the wall normal while :yref:`overlap<VolumeGeom.penetrationVolume>` and :yref:`contact points<VolumeGeom.contactPoint>` are defined likewise to :yref:`Ig2_LevelSet_LevelSet_VolumeGeom`.",
	((uint,nRefineOctree,5,,"The number of refinements performed by the Octree algorithm used to compute the overlap volume between two particles. Default is 5."))
	((Real,smearCoeffOctree,1.0,,"Smearing coefficient for the smeared Heaviside step function in the overlap volume integration. The transition width, or smearing width, is equal to half the diagonal of the smallest integration cell divided by the smearing coefficient."))
	((bool,useAABE,false,,"If true, use the provided (locally) axis-aligned bounding ellipsoid (AABE) to reduce the potential overlap volume between the particles. Increases accuracy of the Octree algrithm because the smallest integration cells will be smaller."))
	,,
	);
	// clang-format on
	DECLARE_LOGGER;
	FUNCTOR2D(Wall, LevelSet);
	DEFINE_FUNCTOR_ORDER_2D(Wall, LevelSet);
};
REGISTER_SERIALIZABLE(Ig2_Wall_LevelSet_VolumeGeom);


class ShopLSvolume {
public:
	static shared_ptr<VolumeGeom> volGeomPtr(
	        Vector3r                       ctctPt,
	        Real                           un,
	        Real                           rad1,
	        Real                           rad2,
	        const State&                   rbp1,
	        const State&                   rbp2,
	        const shared_ptr<Interaction>& c,
	        const Vector3r&                currentNormal,
	        const Vector3r&                shift2);
	static shared_ptr<VolumeGeom> volGeomPtrForLaterRemoval(const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c);
	struct overlapRegionData {
		Real     volume    = 0.0;              // Volume of overlap region
		Vector3r normal1   = Vector3r::Zero(); // Weighted level-set normal particle 1
		Vector3r normal2   = Vector3r::Zero(); // Weighted level-set normal particle 2
		Vector3r centroid1 = Vector3r::Zero(); // Weighted centroid particle 1
		Vector3r centroid2 = Vector3r::Zero(); // Weighted centroid particle 2
		Real     depth1    = 0.0;
		Real     depth2    = 0.0;
		Real     area      = 0.0;
	};
	struct layerData {
		//layerData(Real a = 0.0, Real b = 0.0, Vector3r c = Vector3r::Zero(), Vector3r d = Vector3r::Zero())
		//	: cellVolume(a), Rmax(b), refineStep1(c), refineStep2(d) {}
		Real                  cellVolume = 0.0; // Volume of a single cell at the current refinement level
		Real                  Rmax       = 0.0;
		std::vector<Vector3r> refineStep1 { std::vector<Vector3r>(8, Vector3r::Zero()) };
		std::vector<Vector3r> refineStep2 { std::vector<Vector3r>(8, Vector3r::Zero()) };
	};
	static Vector3r          normalToPointOnEllipsoid(const Vector3r normal, const Vector3r prinAxesSq);
	static overlapRegionData recursiveVolumeIntegration( // Move to ShopLS.*pp or keep here? Maybe nice to have all relevant functions in this file.
	        const shared_ptr<LevelSet>&   lsShape1,
	        const shared_ptr<LevelSet>&   lsShape2,
	        const Vector3r                pt1,
	        const Vector3r                pt2,
	        const std::vector<layerData>& layers, // cell volume, Rmax, refinement steps for grid 1 and 2
	        const uint                    layerId,
	        const Real                    smearCoeffOctree);
	static overlapRegionData recursiveVolumeIntegrationWall(
	        const shared_ptr<LevelSet>&   lsShape,
	        const Vector3r                pt,
	        const Vector3r                nWall,
	        const Real                    kWall,
	        const std::vector<layerData>& layers, // cell volume, Rmax, refinement steps for grid 1 and 2
	        const uint                    layerId,
	        const Real                    smearCoeffOctree);
};

} // namespace yade
#endif // YADE_LS_DEM
