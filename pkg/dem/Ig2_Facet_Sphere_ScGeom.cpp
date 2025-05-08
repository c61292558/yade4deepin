/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Ig2_Facet_Sphere_ScGeom.hpp"
#include <lib/base/Math.hpp>
#include <core/InteractionLoop.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((Ig2_Facet_Sphere_ScGeom)(Ig2_Facet_Sphere_ScGeom6D)(Ig2_Wall_Sphere_ScGeom));

CREATE_LOGGER(Ig2_Facet_Sphere_ScGeom);

bool Ig2_Facet_Sphere_ScGeom::go(
        const shared_ptr<Shape>&       cm1,
        const shared_ptr<Shape>&       cm2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	TIMING_DELTAS_START();
	const Se3r& se31  = state1.se3;
	const Se3r& se32  = state2.se3;
	Facet*      facet = static_cast<Facet*>(cm1.get());
	/* could be written as (needs to be tested):
	 * Vector3r cl=se31.orientation.Conjugate()*(se32.position-se31.position);
	 */
	Matrix3r facetAxisT = se31.orientation.toRotationMatrix();
	Matrix3r facetAxis  = facetAxisT.transpose();
	// local orientation
	Vector3r cl = facetAxis
	        * (se32.position + shift2
	           - se31.position); // "contact line" = branch vector from Facet center (= inscribed circle center) to Sphere center, in facet-local coords

	//
	// BEGIN everything in facet-local coordinates
	//

	Vector3r normal = facet->normal;
	Real     L      = normal.dot(cl);
	if (L < 0) {
		normal = -normal;
		L      = -L;
	}

	Real sphereRadius = static_cast<Sphere*>(cm2.get())->radius;
	if (L > sphereRadius && !c->isReal()
	    && !force) { // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false
		TIMING_DELTAS_CHECKPOINT("Ig2_Facet_Sphere_ScGeom");
		return false;
	}

	Vector3r        cp = cl - L * normal; // in-plane component of branch vector
	const Vector3r* ne = facet->ne;       // edges normals

	Real penetrationDepth = 0;

	Real bm = ne[0].dot(cp); // a trilinear coordinate (https://en.wikipedia.org/wiki/Trilinear_coordinates) of in-plane branch vector
	int  m  = 0;
	for (int i = 1; i < 3; ++i) {
		Real b = ne[i].dot(cp);
		if (bm < b) {
			bm = b;
			m  = i;
		}
	}
	// bm is now the biggest trilinear coordinate of in-plane branch vector
	// index of corresponding edge is m

	Real sh  = sphereRadius * shrinkFactor;
	Real icr = facet->icr - sh;
	if (icr < 0) {
		LOG_WARN("a radius of a facet's inscribed circle less than zero! So, shrinkFactor is too large and would be reduced to zero.");
		shrinkFactor = 0;
		icr          = facet->icr;
		sh           = 0;
	}


	if (bm < icr) // projection of sphere's center is strictly within facet's surface (shrinking considerations left aside)
	// precise proof should come from considering trilinear and barycentric coordinates
	{
		LOG_DEBUG("Contact within facet surface");
		penetrationDepth = sphereRadius - L;
		normal.normalize();
	} else { // projection of sphere's center is outside the facet's surface
		LOG_DEBUG("The 'else' part for Facet-Sphere contact");
		cp = cp + ne[m] * (icr - bm);
		if (cp.dot(ne[(m - 1 < 0) ? 2 : m - 1]) > icr) // contact with vertex m
		                                               //			cp = facet->vertices[m];
			cp = facet->vu[m] * (facet->vl[m] - sh);
		else if (cp.dot(ne[m = (m + 1 > 2) ? 0 : m + 1]) > icr) // contact with vertex m+1
		                                                        //			cp = facet->vertices[(m+1>2)?0:m+1];
			cp = facet->vu[m] * (facet->vl[m] - sh);
		normal    = cl - cp;
		Real norm = normal.norm();
		normal /= norm;
		penetrationDepth = sphereRadius - norm;
	}

	//
	// END everything in facet-local coordinates
	//

	if (penetrationDepth > 0 || c->isReal()) {
		shared_ptr<ScGeom> scm;
		bool               isNew = !c->geom;
		if (c->geom) scm = YADE_PTR_CAST<ScGeom>(c->geom);
		else
			scm = shared_ptr<ScGeom>(new ScGeom());

		normal                = facetAxisT * normal; // in global orientation
		scm->contactPoint     = se32.position + shift2 - (sphereRadius - 0.5 * penetrationDepth) * normal;
		scm->penetrationDepth = penetrationDepth;
		scm->radius1          = (hertzian) ? hertzFac * sphereRadius : 2 * sphereRadius;
		scm->radius2          = sphereRadius;
		if (isNew) c->geom = scm;
		scm->precompute(state1, state2, scene, c, normal, isNew, shift2, false /*avoidGranularRatcheting only for sphere-sphere*/);
		TIMING_DELTAS_CHECKPOINT("Ig2_Facet_Sphere_ScGeom");
		return true;
	}
	TIMING_DELTAS_CHECKPOINT("Ig2_Facet_Sphere_ScGeom");
	return false;
}


bool Ig2_Facet_Sphere_ScGeom::goReverse(
        const shared_ptr<Shape>&       cm1,
        const shared_ptr<Shape>&       cm2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	//LOG_WARN("Swapped interaction order for "<<c->getId2()<<"&"<<c->getId1());
	return go(cm2, cm1, state2, state1, -shift2, force, c);
}

bool Ig2_Facet_Sphere_ScGeom6D::go(
        const shared_ptr<Shape>&       cm1,
        const shared_ptr<Shape>&       cm2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	bool isNew = !c->geom;
	if (Ig2_Facet_Sphere_ScGeom::go(cm1, cm2, state1, state2, shift2, force, c)) {
		if (isNew) { //generate a 6DOF interaction from the 3DOF one generated by Ig2_Facet_Sphere_ScGeom
			shared_ptr<ScGeom6D> sc(new ScGeom6D());
			*(YADE_PTR_CAST<ScGeom>(sc)) = *(YADE_PTR_CAST<ScGeom>(c->geom));
			c->geom                      = sc;
		}
		YADE_PTR_CAST<ScGeom6D>(c->geom)->precomputeRotations(state1, state2, isNew, false);
		return true;
	} else
		return false;
}


bool Ig2_Facet_Sphere_ScGeom6D::goReverse(
        const shared_ptr<Shape>&       cm1,
        const shared_ptr<Shape>&       cm2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	return go(cm2, cm1, state2, state1, -shift2, force, c);
}


/********* Wall + Sphere **********/

bool Ig2_Wall_Sphere_ScGeom::go(
        const shared_ptr<Shape>&       cm1,
        const shared_ptr<Shape>&       cm2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	Wall*      wall   = static_cast<Wall*>(cm1.get());
	const Real radius = static_cast<Sphere*>(cm2.get())->radius;
	const int& ax(wall->axis);
	Real       dist = (state2.pos)[ax] + shift2[ax] - state1.pos[ax];         // signed "distance" between centers
	if (!c->isReal() && math::abs(dist) > radius && !force) { return false; } // wall and sphere too far from each other

	// contact point is sphere center projected onto the wall
	Vector3r contPt = state2.pos + shift2;
	contPt[ax]      = state1.pos[ax];
	Vector3r normal(0., 0., 0.);
	// wall interacting from both sides: normal depends on sphere's position
	assert(wall->sense == -1 || wall->sense == 0 || wall->sense == 1);
	if (wall->sense == 0) normal[ax] = dist > 0 ? 1. : -1.;
	else
		normal[ax] = wall->sense == 1 ? 1. : -1;

	bool isNew = !c->geom;
	if (isNew) c->geom = shared_ptr<ScGeom>(new ScGeom());
	const shared_ptr<ScGeom>& ws = YADE_PTR_CAST<ScGeom>(c->geom);
	ws->radius1          = (hertzian) ? hertzFac * radius : radius; // default (not Hertzian): wall's "radius" is chosen as the same as the sphere's radius
	ws->radius2          = radius;
	ws->contactPoint     = contPt;
	ws->penetrationDepth = -(math::abs(dist) - radius);
	// ws->normal is assigned by precompute
	ws->precompute(state1, state2, scene, c, normal, isNew, shift2, noRatch);
	return true;
}

} // namespace yade
