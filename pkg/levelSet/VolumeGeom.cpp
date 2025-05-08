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
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/levelSet/VolumeGeom.hpp>

namespace yade {
YADE_PLUGIN((VolumeGeom));

void VolumeGeom::precompute(
        const State&                   rbp1,
        const State&                   rbp2,
        const Scene*                   scene,
        const shared_ptr<Interaction>& c,
        const Vector3r&                currentNormal,
        bool                           isNew,
        const Vector3r&                shift2)
{
	if (!isNew) {
		orthonormal_axis = normal.cross(currentNormal);
		Real angle       = scene->dt * 0.5 * normal.dot(rbp1.angVel + rbp2.angVel);
		twist_axis       = angle * normal;
	} else {
		twist_axis = orthonormal_axis = Vector3r::Zero();
	}
	// Update contact normal
	normal = currentNormal;
	// Precompute shear increment
	Vector3r c1x = (contactPoint - rbp1.pos);
	Vector3r c2x;
	if (scene->isPeriodic) { // Ternary operator won't work, would need scoping.
		c2x = contactPoint - (rbp2.pos + shift2);
	} else {
		c2x = contactPoint - rbp2.pos;
	}
	Vector3r relativeVelocity = (rbp2.vel + rbp2.angVel.cross(c2x)) - (rbp1.vel + rbp1.angVel.cross(c1x));
	// Add the shift velocity of the periodic box if periodicity is being used
	if (scene->isPeriodic) { relativeVelocity += scene->cell->intrShiftVel(c->cellDist); }
	// Keep the shear part only
	relativeVelocity = relativeVelocity - normal.dot(relativeVelocity) * normal;
	shearInc         = relativeVelocity * scene->dt;
}


Vector3r& VolumeGeom::rotate(Vector3r& shearForce) const
{
	// Formally, is this correct for non-shperical geometries? There was a rotateNonSpherical in ScGeom.cpp. However, it is not being used. Only in KnKsPBLaw.cpp.
	// Approximated rotations
	shearForce -= shearForce.cross(orthonormal_axis);
	shearForce -= shearForce.cross(twist_axis);
	// NOTE : make sure it is in the tangent plane? It's never been done before. Is it not adding rounding errors at the same time in fact?...
	// shearForce -= normal.dot(shearForce) * normal;
	return shearForce;
}


} // namespace yade
#endif //YADE_LS_DEM
