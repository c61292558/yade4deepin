/*****************************************************************************
*  2023 © DLH van der Haven, dannyvdhaven@gmail.com, University of Cambridge *
*                                                                            *
*  For details, please see van der Haven et al., "A physically consistent    *
*   Discrete Element Method for arbitrary shapes using Volume-interacting    *
*   Level Sets", Comput. Methods Appl. Mech. Engrg., 414 (116165):1-21       *
*   https://doi.org/10.1016/j.cma.2023.116165                                *
*  This project has been financed by Novo Nordisk A/S (Bagsværd, Denmark).   *
*                                                                            *
*  This program is licensed under GNU GPLv2, see file LICENSE for details.   *
*****************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <lib/base/Math.hpp>
#include <core/IGeom.hpp>
#include <core/Interaction.hpp>
#include <core/State.hpp>
#include <pkg/levelSet/LevelSet.hpp>

namespace yade {

class VolumeGeom : public IGeom {
public:
	virtual ~VolumeGeom()         = default;
	VolumeGeom(const VolumeGeom&) = default;
	inline VolumeGeom& operator   =(const VolumeGeom& source)
	{
		penetrationVolume       = source.penetrationVolume;
		contactArea             = source.contactArea;
		averagePenetrationDepth = source.averagePenetrationDepth;
		contactPoint            = source.contactPoint;
		normal                  = source.normal;
		twist_axis              = source.twist_axis;
		orthonormal_axis        = source.orthonormal_axis;
		shearInc                = source.shearInc;
		radius1                 = source.radius1;
		radius2                 = source.radius2;
		return *this;
	}

	// Precompute data for shear evaluation
	void precompute(
	        const State&                   rbp1,
	        const State&                   rbp2,
	        const Scene*                   scene,
	        const shared_ptr<Interaction>& c,
	        const Vector3r&                currentNormal,
	        bool                           isNew,
	        const Vector3r&                shift2);
	Vector3r&       rotate(Vector3r& shearForce) const;
	const Vector3r& shearIncrement() const { return shearInc; };
	Real            refR1, refR2, radius1, radius2;

protected:
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(VolumeGeom,IGeom,
		"Geometry of the interaction between two :yref:`LevelSet` bodies when using volume-based interactions. Will soon become the general class for volume interaction, such that it works for polyhedra as well.",
			((Real,penetrationVolume,NaN,(Attr::noSave|Attr::readonly),"Volume of the overlap or penetrating region."))
			((Real,contactArea,NaN,(Attr::noSave|Attr::readonly),"Contact area perpendicular to the normal."))
			((Real,averagePenetrationDepth,NaN,(Attr::noSave|Attr::readonly),"penetrationVolume / contactArea."))
			((Vector3r,contactPoint,Vector3r::Zero(),,"Contact point (global coordinates), centroid of the penetration volume."))
			((Vector3r,normal,Vector3r::Zero(),,"Normal direction of the interaction."))
			((Vector3r,twist_axis,Vector3r::Zero(),,""))
			((Vector3r,orthonormal_axis,Vector3r::Zero(),,""))
			((Vector3r,shearInc,Vector3r::Zero(),(Attr::noSave|Attr::readonly),"Shear displacement increment in the last step."))
			,
			createIndex();
			twist_axis = Vector3r::Zero();
			orthonormal_axis = Vector3r::Zero();
		);
	// clang-format on
	REGISTER_CLASS_INDEX(VolumeGeom, IGeom);
};
REGISTER_SERIALIZABLE(VolumeGeom);

} // namespace yade
#endif // YADE_LS_DEM
