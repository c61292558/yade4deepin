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
#include <lib/high-precision/Constants.hpp>
#include <pkg/levelSet/LevelSetVolumeIg2.hpp>
#include <pkg/levelSet/ShopLS.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade {
YADE_PLUGIN((Ig2_LevelSet_LevelSet_VolumeGeom)(Ig2_Wall_LevelSet_VolumeGeom)); // To add:(Ig2_Box_LevelSet_VolumeGeom)
//CREATE_LOGGER(Ig2_Box_LevelSet_VolumeGeom); // To be implemented
CREATE_LOGGER(Ig2_LevelSet_LevelSet_VolumeGeom);
CREATE_LOGGER(Ig2_Wall_LevelSet_VolumeGeom);


shared_ptr<VolumeGeom> ShopLSvolume::volGeomPtr(
        Vector3r                       ctctPt,
        Real                           un,
        Real                           rad1,
        Real                           rad2,
        const State&                   rbp1,
        const State&                   rbp2,
        const shared_ptr<Interaction>& c,
        const Vector3r&                currentNormal,
        const Vector3r&                shift2)
{
	shared_ptr<VolumeGeom> volGeomPtr;
	bool                   isNew = !c->geom;
	if (isNew) volGeomPtr = shared_ptr<VolumeGeom>(new VolumeGeom());
	else
		volGeomPtr = YADE_PTR_CAST<VolumeGeom>(c->geom);
	volGeomPtr->contactPoint      = ctctPt;
	volGeomPtr->penetrationVolume = un;
	// NB radius1, radius2: are still useful for time step determination with respect to rotational stiffnesses
	// and also, as refR1, refR2, for contact stiffness expression in FrictPhys/FrictMat.
	volGeomPtr->radius1 = rad1;
	volGeomPtr->radius2 = rad2;
	volGeomPtr->refR1   = rad1;
	volGeomPtr->refR2   = rad2;
	volGeomPtr->precompute(rbp1, rbp2, Omega::instance().getScene().get(), c, currentNormal, isNew, shift2);
	// precompute will take care of
	// * preparing the rotation of shearForce to the new tangent plane (done later, in Law2) defining these orthonormal_axis and twist_axis
	// * updating geomPtr->normal to normal
	// * computing the relative velocity at contact, through getIncidentVel(avoidGranularRatcheting=false), using now-defined contactPoint
	return volGeomPtr;
}

shared_ptr<VolumeGeom> ShopLSvolume::volGeomPtrForLaterRemoval(const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c)
{
	// Use when no overlap, but still need to have some geom data (while returning true for InteractionLoop workflow).
	// Otherwise we would need to update InteractionLoop itself to avoid LOG_WARN messages).
	// Data mostly include an infinite tensile stretch to insure subsequent interaction removal (by Law2).
	return ShopLSvolume::volGeomPtr(
	        Vector3r::Zero() /* Inconsequential bullsh..*/,
	        -std::numeric_limits<Real>::infinity() /* Arbitrary big tensile value to trigger interaction removal by Law2*/,
	        1, /* Inconsequential bullsh..*/
	        1, /* Inconsequential bullsh..*/
	        rbp1,
	        rbp2,
	        c,
	        Vector3r::UnitX() /* Inconsequential bullsh..*/,
	        Vector3r::Zero() /* Inconsequential bullsh..*/);
}


/*************************
*   PARTICLE-PARTICLE    *
**************************/

bool Ig2_LevelSet_LevelSet_VolumeGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	// Determine the AABB zone where the bodies' bounds overlap. TODO: possible use of Eigen AlignedBox?
	std::array<Real, 6>     overlapAABB; // Format: [xmin,xmax,ymin,ymax,zmin,zmax]
	const shared_ptr<Bound> bound1 = Body::byId(c->id1, scene)->bound;
	const shared_ptr<Bound> bound2 = Body::byId(c->id2, scene)->bound;
	for (unsigned int axis = 0; axis < 3; axis++) {
		overlapAABB[2 * axis]     = math::max(bound1->min[axis], bound2->min[axis] + shift2[axis]);
		overlapAABB[2 * axis + 1] = math::min(bound1->max[axis], bound2->max[axis] + shift2[axis]);
		if (overlapAABB[2 * axis + 1]
		    < overlapAABB[2 * axis]) { // Overlap AABB is empty. Possible when bodies' bounds themselves are parallel and overlapping.
			if (c->isReal()) {     // Check if the potential interaction is new in this step.
				// Mark current interaction for removal.
				c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
				return true;
			} else {
				// Do not create new interaction.
				return false;
			}
		}
	}
	Vector3r minBoOverlap  = Vector3r(overlapAABB[0], overlapAABB[2], overlapAABB[4]),
	         maxBoOverlap  = Vector3r(overlapAABB[1], overlapAABB[3], overlapAABB[5]);
	Vector3r sizeBoOverlap = maxBoOverlap - minBoOverlap;
	Real     volBoOverlap  = sizeBoOverlap.prod();
	Vector3r cntrBoOverlap = minBoOverlap + sizeBoOverlap / 2;
	// The AABB size and centre can be randomised a bit with (2*double(rand())/RAND_MAX-1)*sizeBo/pow(8,nRefineOctree-1)*4
	// in attempt to smear numerical error, but the effect is not big.

	// Positions and orientations of the particles at start and end
	shared_ptr<LevelSet> lsShape1 = YADE_PTR_CAST<LevelSet>(shape1);
	shared_ptr<LevelSet> lsShape2 = YADE_PTR_CAST<LevelSet>(shape2);
	Vector3r             centrEnd1(state1.pos), centrEnd2(state2.pos + shift2);
	// ori = rotation from reference configuration (local axes) to current one (global axes)
	// ori.conjugate() from the current configuration (global axes) to the reference one (local axes)
	Quaternionr rot1(state1.ori), rot2(state2.ori), rotConj1(state1.ori.conjugate()), rotConj2(state2.ori.conjugate());

	// Box reduction based on projection on ellipsoid slices.
	Quaternionr rotGuess1, rotGuess2, rotGuess1conj, rotGuess2conj;
	if (useAABE && lsShape1->hasAABE && lsShape2->hasAABE) {
		// These properties should be added to the level set shape and retrieved here with something like
		// Vector3r prinAxes1 = lsShape1->principalAxes();
		Vector3r prinAxes1 = lsShape1->axesAABE, prinAxes2 = lsShape2->axesAABE;

		// Check if the bounding spheres overlap. If not, do not create an interaction.
		Vector3r centrToCentr = centrEnd2 - centrEnd1;
		Real     ctcDist      = centrToCentr.norm();
		if (!c->isReal() && !force && (ctcDist > (prinAxes1.maxCoeff() + prinAxes2.maxCoeff()))) {
			return false;
		} else if (ctcDist > (prinAxes1.maxCoeff() + prinAxes2.maxCoeff())) {
			// Mark current interaction for removal.
			c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
			return true;
		}

		// Get a sensible guess for the normal of the interaction between the particles.
		Vector3r nGuess      = centrToCentr.normalized();                          // Normal, global coordinates, from 1 to 2
		Vector3r nzGuessLoc1 = rotConj1 * nGuess, nzGuessLoc2 = rotConj2 * nGuess; // Normal, local coordinates, from 1 to 2

		// Squared principal half axes
		Vector3r axSq1 = Vector3r(prinAxes1[0] * prinAxes1[0], prinAxes1[1] * prinAxes1[1], prinAxes1[2] * prinAxes1[2]),
		         axSq2 = Vector3r(prinAxes2[0] * prinAxes2[0], prinAxes2[1] * prinAxes2[1], prinAxes2[2] * prinAxes2[2]);

		// To bound the potential overlap, find the farthest point on the ellipsoid in the direction of the guessed interaction normal
		Vector3r farPtz1 = ShopLSvolume::normalToPointOnEllipsoid(nzGuessLoc1, axSq1);
		Vector3r farPtz2 = -ShopLSvolume::normalToPointOnEllipsoid(nzGuessLoc2, axSq2);

		// We now proceed to make a bounding box for the overlap volume between the two particles,
		// based on the axis-aligned bounding ellipsoids. To do this, we will rotate the ellipsoids
		// in local coordinates to align with the guessed interation normal. This allows us to
		// easily slice the ellipsoids at the farthest point in the direction of the guessed normal.

		// Get the rotation that moves the local guessed interaction normal onto the local z axis (and back)
		rotGuess1     = Quaternionr().setFromTwoVectors(nzGuessLoc1, Vector3r::UnitZ()),
		rotGuess2     = Quaternionr().setFromTwoVectors(nzGuessLoc2, Vector3r::UnitZ());
		rotGuess1conj = rotGuess1.conjugate(), rotGuess2conj = rotGuess2.conjugate();

		// Rotate farthest points, aling interaction normal onto the z axis
		Vector3r farPtz1rot = rotGuess1 * farPtz1, farPtz2rot = rotGuess2 * farPtz2;

		// Check if the ellipsoids can touch by using the farthest z coordinates of the ellipsoids
		if ((farPtz2rot[2] + ctcDist - farPtz1rot[2]) < 0.0) {
			// The matrix describing a general ellipsoid (not restricted to the axis-aligned case)
			Matrix3r A1 = Vector3r(1.0 / axSq1[0], 1.0 / axSq1[1], 1.0 / axSq1[2]).asDiagonal(),
			         A2 = Vector3r(1.0 / axSq2[0], 1.0 / axSq2[1], 1.0 / axSq2[2]).asDiagonal();

			// Rotate the ellipsoids using the tensor relation: Anew = R*A*inv(R)
			// The coordinate systems of the ellipsoids are now aliged but with a z offset of ctcDist
			Matrix3r A1rot = rotGuess1 * A1 * rotGuess1.inverse(), A2rot = rotGuess2 * A2 * rotGuess2.inverse();

			// We now do a plane-ellipsoid intersection at the extreme coordinates as given by the
			// farthest points on the ellipsoid in direction of the guessed interaction normal.
			// The resulting extremes of the resulting ellipses will give us the x and y limits of
			// the bounding box for the potential interaction.

			// Should funtionalise this at some point.
			// The full equation of the ellipse resulting from a sliced ellipsoid at z = h:
			// 	axx x^2 + ayy y^2 + azz h^2 + 2 axy x y + 2 axz x h + 2 ayz y h = 1
			// The h are the extreme z coordinates and locations to intersect the ellipsoids at
			Real axx = A1rot(0, 0), ayy = A1rot(1, 1), azz = A1rot(2, 2), axy = A1rot(0, 1), axz = A1rot(0, 2), ayz = A1rot(1, 2);

			// Location on z axis for the slicing plane, determined by extreme of the OTHER ellipsoid
			// (because that is the one that penetrates).
			// NB: For the rotated ellipsoids, particle 1 is always below particle 2!!
			Real h = farPtz2rot[2] + ctcDist;

			// Compute x limits of ellipse
			Real xterm1 = 0.5 / (axx - axy * axy / ayy), xterm2 = 2 * axz * h - (2 * axy * ayz * h) / ayy,
			     xterm3 = sqrt(pow(-xterm2, 2) - 4 * (-axy * axy / ayy + axx) * (-1 - ayz * ayz * h * h / ayy + azz * h * h));
			Real xlim11 = xterm1 * (xterm2 - xterm3), xlim12 = xterm1 * (xterm2 + xterm3);
			Real xlim1neg, xlim1pos;
			if (xlim11 > xlim12) {
				xlim1pos = xlim11;
				xlim1neg = xlim12;
			} else {
				xlim1pos = xlim12;
				xlim1neg = xlim11;
			}

			// Compute y limits of ellipse
			Real yterm1 = 0.5 / (ayy - axy * axy / axx), yterm2 = 2 * ayz * h - (2 * axy * axz * h) / axx,
			     yterm3 = sqrt(pow(-yterm2, 2) - 4 * (-axy * axy / axx + ayy) * (-1 - axz * axz * h * h / axx + azz * h * h));
			Real ylim11 = yterm1 * (yterm2 + yterm3), ylim12 = yterm1 * (yterm2 - yterm3);
			Real ylim1neg, ylim1pos;
			if (ylim11 > ylim12) {
				ylim1pos = ylim11;
				ylim1neg = ylim12;
			} else {
				ylim1pos = ylim12;
				ylim1neg = ylim11;
			}

			// Check if the slice we took is not actually smaller than the widest slice
			// within the potential overlap volume. Check along x axis.
			Vector3r nxGuessLoc1 = rotGuess1conj * Vector3r::UnitX();
			Vector3r farPtx1     = ShopLSvolume::normalToPointOnEllipsoid(nxGuessLoc1, axSq1);
			Vector3r farPtx1rot  = rotGuess1 * farPtx1;
			if (farPtx1rot[2] > h) { // Point is above slice, needs to be considered
				if (farPtx1rot[0] > 0) {
					xlim1pos = std::max(xlim1pos, farPtx1rot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, farPtx1rot[0]);
				}
				if (farPtx1rot[1] > 0) {
					ylim1pos = std::max(ylim1pos, farPtx1rot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, farPtx1rot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			// is -farPt1xrot due to inversion symmetry.
			if (-farPtx1rot[2] > h) { // Point is above slice, needs to be considered
				if (-farPtx1rot[0] > 0) {
					xlim1pos = std::max(xlim1pos, -farPtx1rot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, -farPtx1rot[0]);
				}
				if (-farPtx1rot[1] > 0) {
					ylim1pos = std::max(ylim1pos, -farPtx1rot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, -farPtx1rot[1]);
				}
			}
			// Now check along y axis.
			Vector3r nyGuessLoc1 = rotGuess1conj * Vector3r::UnitY();
			Vector3r farPty1     = ShopLSvolume::normalToPointOnEllipsoid(nyGuessLoc1, axSq1);
			Vector3r farPty1rot  = rotGuess1 * farPty1;
			if (farPty1rot[2] > h) { // Point is above slice, needs to be considered
				if (farPty1rot[0] > 0) {
					xlim1pos = std::max(xlim1pos, farPty1rot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, farPty1rot[0]);
				}
				if (farPty1rot[1] > 0) {
					ylim1pos = std::max(ylim1pos, farPty1rot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, farPty1rot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			// is -farPt1xrot due to inversion symmetry.
			if (-farPty1rot[2] > h) { // Point is above slice, needs to be considered
				if (-farPty1rot[0] > 0) {
					xlim1pos = std::max(xlim1pos, -farPty1rot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, -farPty1rot[0]);
				}
				if (-farPty1rot[1] > 0) {
					ylim1pos = std::max(ylim1pos, -farPty1rot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, -farPty1rot[1]);
				}
			}

			//
			// Repeat for the other particle
			//
			axx = A2rot(0, 0), ayy = A2rot(1, 1), azz = A2rot(2, 2), axy = A2rot(0, 1), axz = A2rot(0, 2), ayz = A2rot(1, 2);
			h = farPtz1rot[2] - ctcDist; // Other way around for particle 2

			// Compute x limits of ellipse
			xterm1 = 0.5 / (axx - axy * axy / ayy), xterm2 = 2 * axz * h - (2 * axy * ayz * h) / ayy,
			xterm3      = sqrt(pow(-xterm2, 2) - 4 * (-axy * axy / ayy + axx) * (-1 - ayz * ayz * h * h / ayy + azz * h * h));
			Real xlim21 = xterm1 * (xterm2 - xterm3), xlim22 = xterm1 * (xterm2 + xterm3);
			Real xlim2neg, xlim2pos;
			if (xlim21 > xlim22) {
				xlim2pos = xlim21;
				xlim2neg = xlim22;
			} else {
				xlim2pos = xlim22;
				xlim2neg = xlim21;
			}

			// Compute y limits of ellipse
			yterm1 = 0.5 / (ayy - axy * axy / axx), yterm2 = 2 * ayz * h - (2 * axy * axz * h) / axx,
			yterm3      = sqrt(pow(-yterm2, 2) - 4 * (-axy * axy / axx + ayy) * (-1 - axz * axz * h * h / axx + azz * h * h));
			Real ylim21 = yterm1 * (yterm2 + yterm3), ylim22 = yterm1 * (yterm2 - yterm3);
			Real ylim2neg, ylim2pos;
			if (ylim21 > ylim22) {
				ylim2pos = ylim21;
				ylim2neg = ylim22;
			} else {
				ylim2pos = ylim22;
				ylim2neg = ylim21;
			}

			// Check along x axis.
			Vector3r nxGuessLoc2 = rotGuess2conj * Vector3r::UnitX();
			Vector3r farPtx2     = ShopLSvolume::normalToPointOnEllipsoid(nxGuessLoc2, axSq2);
			Vector3r farPtx2rot  = rotGuess2 * farPtx2;
			if (farPtx2rot[2] < h) { // Point is below slice, needs to be considered
				if (farPtx2rot[0] > 0) {
					xlim2pos = std::max(xlim2pos, farPtx2rot[0]);
				} else {
					xlim2neg = std::min(xlim2neg, farPtx2rot[0]);
				}
				if (farPtx2rot[1] > 0) {
					ylim2pos = std::max(ylim2pos, farPtx2rot[1]);
				} else {
					ylim2neg = std::min(ylim2neg, farPtx2rot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			if (-farPtx2rot[2] < h) { // Point is below slice, needs to be considered
				if (-farPtx2rot[0] > 0) {
					xlim2pos = std::max(xlim2pos, -farPtx2rot[0]);
				} else {
					xlim2neg = std::min(xlim2neg, -farPtx2rot[0]);
				}
				if (-farPtx2rot[1] > 0) {
					ylim2pos = std::max(ylim2pos, -farPtx2rot[1]);
				} else {
					ylim2neg = std::min(ylim2neg, -farPtx2rot[1]);
				}
			}
			// Now check along y axis.
			Vector3r nyGuessLoc2 = rotGuess2conj * Vector3r::UnitY();
			Vector3r farPty2     = ShopLSvolume::normalToPointOnEllipsoid(nyGuessLoc2, axSq2);
			Vector3r farPty2rot  = rotGuess2 * farPty2;
			if (farPty2rot[2] < h) { // Point is below slice, needs to be considered
				if (farPty2rot[0] > 0) {
					xlim2pos = std::max(xlim2pos, farPty2rot[0]);
				} else {
					xlim2neg = std::min(xlim2neg, farPty2rot[0]);
				}
				if (farPty2rot[1] > 0) {
					ylim2pos = std::max(ylim2pos, farPty2rot[1]);
				} else {
					ylim2neg = std::min(ylim2neg, farPty2rot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			if (-farPty2rot[2] < h) { // Point is below slice, needs to be considered
				if (-farPty2rot[0] > 0) {
					xlim2pos = std::max(xlim2pos, -farPty2rot[0]);
				} else {
					xlim2neg = std::min(xlim2neg, -farPty2rot[0]);
				}
				if (-farPty2rot[1] > 0) {
					ylim2pos = std::max(ylim2pos, -farPty2rot[1]);
				} else {
					ylim2neg = std::min(ylim2neg, -farPty2rot[1]);
				}
			}

			// Bounding box dimensions in local rotated frame of particle 1, intersecting the limits set by each ellipsoid
			Vector3r minBoOverlap1 = Vector3r(std::max(xlim1neg, xlim2neg), std::max(ylim1neg, ylim2neg), farPtz2rot[2] + ctcDist);
			Vector3r maxBoOverlap1 = Vector3r(std::min(xlim1pos, xlim2pos), std::min(ylim1pos, ylim2pos), farPtz1rot[2]);

			// Initialise limits in global coordinates
			Vector3r minBoOverlapOld = Vector3r(
			        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
			Vector3r maxBoOverlapOld = Vector3r(
			        -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
			// Move the old global AABB to local coordinates
			for (int i = 0; i < 2; i++) { // Go through all possible corners
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						// Select one of the corners of the AABB
						Vector3r cornerGlob = Vector3r(
						        minBoOverlap[0] + i * (sizeBoOverlap[0]),
						        minBoOverlap[1] + j * (sizeBoOverlap[1]),
						        minBoOverlap[2] + k * (sizeBoOverlap[2]));
						// Move corner to local coordinates
						Vector3r cornerLoc = rotGuess1 * rotConj1 * (cornerGlob - centrEnd1);
						// Find maximum and minimum values of the corners
						for (int m = 0; m < 3; m++) {
							if (cornerLoc[m] > maxBoOverlapOld[m]) { maxBoOverlapOld[m] = cornerLoc[m]; }
							if (cornerLoc[m] < minBoOverlapOld[m]) { minBoOverlapOld[m] = cornerLoc[m]; }
						}
					}
				}
			}

			// Recompute the AABB by intersecting the individual bounding boxes
			// All still local and z-aligned!
			for (int i = 0; i < 3; i++) {
				minBoOverlap1[i] = math::max(minBoOverlap1[i], minBoOverlapOld[i]);
				maxBoOverlap1[i] = math::min(maxBoOverlap1[i], maxBoOverlapOld[i]);
			}

			// Compute the final bounding box (BB) parameters.
			// Note: Depending on how tight the ellipsoid fits around the shape,
			// the BB can be a tight fit as well. This means that a litte bit of volume that would
			// be added because of smearing is not added. Currently, it seems that this works best.
			Vector3r cntrBoOverlap1 = (minBoOverlap1 + maxBoOverlap1) / 2;
			cntrBoOverlap           = rot1 * (rotGuess1conj * cntrBoOverlap1) + centrEnd1;
			sizeBoOverlap           = (maxBoOverlap1 - minBoOverlap1);
			volBoOverlap            = sizeBoOverlap.prod();
		} else {
			volBoOverlap = 0.0;
		}

		// Test if there is any overlap at all. Otherwise, do not create an interaction.
		if (!c->isReal() && !force && volBoOverlap < Mathr::EPSILON) {
			return false;
		} else if (volBoOverlap < Mathr::EPSILON) {
			// Mark current interaction for removal.
			c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
			return true;
		}
	}
	// End of overlap bounding box computations

	// Determine the initial query point (center of the overlap bounding box) for the volume intergration in local axes
	// Doing this avoids that we have to do a rigid mapping for every new query point upon refining the intergration mesh
	// The disadvantage is that we have to compute the location of each point twice, but this is still cheaper than ShopLS::rigidMapping.
	Vector3r pt1 = rotConj1 * (cntrBoOverlap - centrEnd1), pt2 = rotConj2 * (cntrBoOverlap - centrEnd2);

	// Precompute a number of parameters that depend on the refinement (number of layers) of the integration mesh.
	// Start at the most coarse layer and end at the finest layer. Note that the refinement step has to be rotated as the
	// orientation of the overlap AABB is different in the respective local axes.
	std::vector<ShopLSvolume::layerData> layers(nRefineOctree); //was layers
	Real RmaxBoOverlap = 0.5 * sqrt(sizeBoOverlap[0] * sizeBoOverlap[0] + sizeBoOverlap[1] * sizeBoOverlap[1] + sizeBoOverlap[2] * sizeBoOverlap[2]);

	if (useAABE) {
		for (uint h = 0; h < nRefineOctree; h++) {
			Real     a           = 1.0 / pow(2.0, h); // Is pow() the correct function to use or is there a high-precision variant?
			Vector3r refineStep  = 0.25 * a * sizeBoOverlap;
			layers[h].cellVolume = volBoOverlap / pow(8, h); // Cell volume
			layers[h].Rmax       = a * RmaxBoOverlap;        // Maximum center-to-edge distance
			int q                = 0;                        // Looping through all the octants
			for (int i = -1; i < 2; i = i + 2) {             // Ugly way to go through all possible sign combinations
				for (int j = -1; j < 2; j = j + 2) {
					for (int k = -1; k < 2; k = k + 2) {
						Vector3r refStepPerm = Vector3r(i * refineStep[0], j * refineStep[1], k * refineStep[2]);
						// Need a different rotation because the AABB is in local coordinates of particle 1
						layers[h].refineStep1[q] = rotGuess1conj * refStepPerm; // Mesh refinement step grid 1
						layers[h].refineStep2[q] = rotGuess2conj * refStepPerm; // Mesh refinement step grid 2
						q++;
					}
				}
			}
		}
	} else {
		for (uint h = 0; h < nRefineOctree; h++) {
			Real     a           = 1.0 / pow(2.0, h); // Is pow() the correct function to use or is there a high-precision variant?
			Vector3r refineStep  = 0.25 * a * sizeBoOverlap;
			layers[h].cellVolume = volBoOverlap / pow(8, h); // Cell volume
			layers[h].Rmax       = a * RmaxBoOverlap;        // Maximum center-to-edge distance
			int q                = 0;                        // Looping through all the octants
			for (int i = -1; i < 2; i = i + 2) {             // Ugly way to go through all possible sign combinations
				for (int j = -1; j < 2; j = j + 2) {
					for (int k = -1; k < 2; k = k + 2) {
						Vector3r refStepPerm     = Vector3r(i * refineStep[0], j * refineStep[1], k * refineStep[2]) + cntrBoOverlap;
						layers[h].refineStep1[q] = rotConj1 * (refStepPerm - centrEnd1) - pt1; // Mesh refinement step grid 1
						layers[h].refineStep2[q] = rotConj2 * (refStepPerm - centrEnd2) - pt2; // Mesh refinement step grid 2
						q++;
					}
				}
			}
		}
	}

	// Recursive volume integration
	// NB: normals and centroids are not yet normalised by the total overlap volume. Also, normals and centroids are still in local coordinates.
	uint                            layerId = 0; // Starting layer, always zero. Needed for recursion.
	ShopLSvolume::overlapRegionData overlap = ShopLSvolume::recursiveVolumeIntegration(lsShape1, lsShape2, pt1, pt2, layers, layerId, smearCoeffOctree);
	Real                            vn      = overlap.volume;

	// Test if there is any overlap at all. Otherwise, do not create an interaction.
	if (!c->isReal() && !force && vn < Mathr::EPSILON) { // Smallest output overlap volume will be just above volAcc/8.
		return false;
	} else if (vn <= Mathr::EPSILON) {
		// Mark current interaction for removal.
		c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
		return true;
	}

	// Centroids and normal vector pointing from particle 1 to particle 2.
	// Normals need not be the same for the particles so we average them.
	Vector3r normal = ((rot1 * overlap.normal1).normalized() - (rot2 * overlap.normal2).normalized()).normalized(), centroidLocal1 = overlap.centroid1 / vn,
	         centroidLocal2 = overlap.centroid2 / vn;
	// Moving centroids to global coordinates
	Vector3r centroid1 = rot1 * centroidLocal1 + centrEnd1, centroid2 = rot2 * centroidLocal2 + centrEnd2;
	// The contact point is the centroid (or center of mass) of the overlap volume
	// The two centers should be exactly the same but we average them just to be sure
	// and reduce numerical error.
	Vector3r contactPoint = 0.5 * (centroid1 + centroid2);

	// Compute a crude estimate of the overlap distance
	Real overlapDistance = pow(vn, 1.0 / 3.0);

	// Compute the center-to-contact distance for the torque computation.
	// NB: Only works if the normal aligns with the center-to-contact vector. We add overlapDistance/2
	// because it is later subtracted when using e.g. ElasticContactLaw (see Law2_VolumeGeom_FrictPhys_Elastic::go in .cpp file).
	// This variable is no longer used because the torque computation is fixed to use the proper center-to-contact vector.
	Real radius1 = (contactPoint - centrEnd1).norm() + 0.5 * overlapDistance, radius2 = (contactPoint - centrEnd2).norm() + 0.5 * overlapDistance;

	// Finalise computation and create or update the interaction
	c->geom = ShopLSvolume::volGeomPtr(contactPoint, vn, radius1, radius2, state1, state2, c, normal, shift2);
	return true;
}

Vector3r ShopLSvolume::normalToPointOnEllipsoid(const Vector3r normal, const Vector3r prinAxesSq)
{
	// Given a normal and an axis-aligned ellipsoid, compute the point on the ellipsoid that has the same normal.
	// This point also equals the farthest point on the ellipsoid in the direction of the normal.
	// Vector prinAxesSq holds the squared principal half axes of the ellipsoid, i.e. (Rx^2,Ry^2,Rz^2).

	/* 
	The solutions are found using (1) the ellipsoid equation (x/Rx)^2 + (y/Ry)^2 + (z/Rz)^2 = 1 
	and (2) the ellipsoid normal n = (2x/Rx^2,2y/Ry^2,2z/Rz^2)/|(2x/Rx^2,2y/Ry^2,2z/Rz^2)|. 
	Solve (1) for x. Substitute found x in ny = 2y/Ry^2/|(2x/Rx^2,2y/Ry^2,2z/Rz^2)|.
	Substitute found x and y in nz = 2z/Rz^2/|(2x/Rx^2,2y/Ry^2,2z/Rz^2)| to get z. 
	Similar procedure for the y component. The x component is then given by the very first
	expression we found for x.
	*/

	// Precompute common term
	Real term1 = prinAxesSq[0] + (prinAxesSq[1] - prinAxesSq[0]) * normal[1] * normal[1] + (prinAxesSq[2] - prinAxesSq[0]) * normal[2] * normal[2];

	// Get y and z components, set x component to zero for now.
	Vector3r pt = Vector3r(
	        0.0,
	        math::sign(normal[1]) * sqrt(normal[1] * normal[1] * prinAxesSq[1] * prinAxesSq[1] / term1),
	        math::sign(normal[2]) * sqrt(normal[2] * normal[2] * prinAxesSq[2] * prinAxesSq[2] / term1));

	// Precompute term
	Real term2 = prinAxesSq.prod() - prinAxesSq[0] * prinAxesSq[2] * pt[1] * pt[1] - prinAxesSq[0] * prinAxesSq[1] * pt[2] * pt[2];

	// Check if the computed term is close to zero, if so, x is also zero. Set x component.
	if (term2 > Mathr::EPSILON) {
		pt[0] = math::sign(normal[0]) * sqrt(term2 / (prinAxesSq[1] * prinAxesSq[2]));
	} else {
		pt[0] = 0.0;
	}

	// Give point on ellipsoid as output.
	// Due to inversion symmetry of axis-aligned ellipsoids, the point with the opposite normal is simply -pt.
	return pt;
}

ShopLSvolume::overlapRegionData
ShopLSvolume::recursiveVolumeIntegration( // Move to ShopLS.*pp or keep here? Maybe nice to have all relevant functions in this file.
        const shared_ptr<LevelSet>&   lsShape1,
        const shared_ptr<LevelSet>&   lsShape2,
        const Vector3r                pt1,
        const Vector3r                pt2,
        const std::vector<layerData>& layers, // cell volume, Rmax, refinement steps for grid 1 and 2
        const uint                    layerId,
        const Real                    smearCoeffOctree)
{
	// A recursive or hierarchical method for integrating the overlap volume of two intersecting level set shapes.
	// The approach is similar to the Octree data structure or adaptive mesh refinement (AMR).
	// NB: Because we allow for concave geometries, the level-set normal at the contact point need not be the correct normal
	// for computing the forces, we thus compute the normal as the volume-weighted normals of the overlap volume segments.

	// Evaluate level set values at query point
	const Real lev1 = lsShape1->distance(pt1, true), lev2 = lsShape2->distance(pt2, true);

	// Structure for output data
	// Mabye working more with pointers would improve the code performance?
	overlapRegionData overlap;

	// Get the details of the current refinement layer
	layerData curLayer = layers[layerId];
	Real      Rmax     = curLayer.Rmax;

	// NB: lev < 0 means we are inside the particle.
	if ((lev1 > Rmax) || (lev2 > Rmax)) {
		// END MESH REFINEMENT
		// Entire region not part of the overlap volume
		return overlap;
	} else if ((lev1 < -Rmax) && (lev2 < -Rmax)) {
		// END MESH REFINEMENT
		// Entire region is part of the overlap volume, add entire cell volume to overlap volume.
		// No normal to add because there is no surface.
		Real volume       = curLayer.cellVolume; // Not sure if copying the volume in a new variable is actually slower or faster.
		overlap.volume    = volume;
		overlap.centroid1 = volume * pt1;
		overlap.centroid2 = volume * pt2;
		overlap.depth1    = volume * lev1;
		overlap.depth2    = volume * lev2;
		return overlap;
	} else {
		// Region has the potential to partly belong to the overlap volume
		if (layerId < layers.size() - 1) { // I find using layers.size() ugly but I stuggled finding a more elegant end condition for the recursion.
			// REFINE THE MESH
			// Generate a new query points, dividing the region in 8 parts of equal size.
			for (int i = 0; i < 8; i++) {
				Vector3r          newPt1 = pt1 + curLayer.refineStep1[i];
				Vector3r          newPt2 = pt2 + curLayer.refineStep2[i];
				overlapRegionData overlapSub
				        = recursiveVolumeIntegration(lsShape1, lsShape2, newPt1, newPt2, layers, layerId + 1, smearCoeffOctree);
				overlap.volume += overlapSub.volume;
				overlap.normal1 += overlapSub.normal1;
				overlap.normal2 += overlapSub.normal2;
				overlap.centroid1 += overlapSub.centroid1;
				overlap.centroid2 += overlapSub.centroid2;
				overlap.depth1 += overlapSub.depth1;
				overlap.depth2 += overlapSub.depth2;
				overlap.area += overlapSub.area;
			}
			return overlap;
		} else {
			// END MESH REFINEMENT
			// Maximum number of refinements reached, only add volume of the cell is the centre is within both particles.
			// This assumption is reasonable because for randomly planes intersecting a rectangle, the expected volume
			// fraction of the box that is overlap is >74% if both centres as inside the box and <14% otherwise.
			// Setting epsilon larger than Rmax leads to errors because larger cells would also start
			// to experience smearing, but smearing is only applied on the smallest cells.
			Real phi = std::max(lev1, lev2);
			Real phiRef(0.);
			if (smearCoeffOctree != 0) phiRef = Rmax / smearCoeffOctree; //Mathr::EPSILON
			if (abs(phi) < phiRef) {
				// Compute the smeared Heaviside step function
				Real volume       = 0.5 * (1.0 - phi / phiRef + sin(-M_PI * phi / phiRef) / M_PI) * curLayer.cellVolume;
				overlap.volume    = volume;
				overlap.centroid1 = volume * pt1;
				overlap.centroid2 = volume * pt2;
				overlap.depth1    = volume * lev1;
				overlap.depth2    = volume * lev2;
				Real surface      = 0;
				Real dS           = pow(curLayer.cellVolume, 2.0 / 3.0);
				if (abs(lev1) < phiRef) {
					Real dS1        = 0.5 / phiRef * (1 + cos(M_PI * lev1 / phiRef)) * dS;
					overlap.normal1 = dS1 * (lsShape1->normal(pt1, true));
					surface         = surface + dS1;
				}
				if (abs(lev2) < phiRef) {
					Real dS2        = 0.5 / phiRef * (1 + cos(M_PI * lev2 / phiRef)) * dS;
					overlap.normal2 = dS2 * (lsShape2->normal(pt2, true));
				}
				overlap.area = surface;
				return overlap;
			} else if (phi < 0) {
				// No normal to add because there is no surface.
				Real volume       = curLayer.cellVolume;
				overlap.volume    = volume;
				overlap.centroid1 = volume * pt1;
				overlap.centroid2 = volume * pt2;
				overlap.depth1    = volume * lev1;
				overlap.depth2    = volume * lev2;
				return overlap;
			} else {
				return overlap;
			}
		}
	} // End mesh refinement check
}


/*************************
*     WALL-PARTICLE      *
**************************/

bool Ig2_Wall_LevelSet_VolumeGeom::go(
        const shared_ptr<Shape>&       shape1,
        const shared_ptr<Shape>&       shape2,
        const State&                   state1,
        const State&                   state2,
        const Vector3r&                shift2,
        const bool&                    force,
        const shared_ptr<Interaction>& c)
{
	// NB: Can and will go wrong if the particle centre penetrates the wall.

	// Particle and wall shapes and positions along axes perpendicular to the wall
	shared_ptr<Wall>     shapeWall = YADE_PTR_CAST<Wall>(shape1);
	shared_ptr<LevelSet> shapeLS   = YADE_PTR_CAST<LevelSet>(shape2);
	Vector3r             wallPos = state1.pos, lsPos = state2.pos + shift2;
	Real                 lsPosAxis(lsPos[shapeWall->axis]), wallPosAxis(wallPos[shapeWall->axis]);

	// Determine the AABB zone where the bodies' bounds overlap. TODO: possible use of Eigen AlignedBox?
	std::array<Real, 6>     overlapAABB; // Format: [xmin,xmax,ymin,ymax,zmin,zmax]
	const shared_ptr<Bound> bound2 = Body::byId(c->id2, scene)->bound;
	for (int axis = 0; axis < 3; axis++) {
		if (axis == shapeWall->axis) {
			if (wallPosAxis >= lsPosAxis) {
				overlapAABB[2 * axis]     = wallPosAxis;
				overlapAABB[2 * axis + 1] = bound2->max[axis] + shift2[axis];
			} else {
				overlapAABB[2 * axis]     = bound2->min[axis] + shift2[axis];
				overlapAABB[2 * axis + 1] = wallPosAxis;
			}
		} else {
			overlapAABB[2 * axis]     = bound2->min[axis] + shift2[axis];
			overlapAABB[2 * axis + 1] = bound2->max[axis] + shift2[axis];
		}

		if (overlapAABB[2 * axis + 1]
		    < overlapAABB[2 * axis]) { // Overlap AABB is empty. Possible when bodies' bounds themselves are parallel and overlapping.
			if (c->isReal()) {     // Check if the potential interaction is new in this step.
				// Mark current interaction for removal.
				c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
				return true;
			} else {
				// Do not create new interaction.
				return false;
			}
		}
	}
	Vector3r minBoOverlap  = Vector3r(overlapAABB[0], overlapAABB[2], overlapAABB[4]),
	         maxBoOverlap  = Vector3r(overlapAABB[1], overlapAABB[3], overlapAABB[5]);
	Vector3r sizeBoOverlap = maxBoOverlap - minBoOverlap;
	Real     volBoOverlap  = sizeBoOverlap.prod();
	Vector3r cntrBoOverlap = minBoOverlap + sizeBoOverlap / 2;
	// End of overlap AABB computations

	// Positions and orientations of the particle and wall at start and end
	Vector3r centrEndLS(lsPos);
	// ori = rotation from reference configuration (local axes) to current one (global axes)
	// ori.conjugate() from the current configuration (global axes) to the reference one (local axes)
	Quaternionr rotLS(state2.ori), rotConjLS(state2.ori.conjugate());

	// Wall normal, NB: Only orientations that align with one of the global axes are allowed.
	Vector3r wallNormal(Vector3r::Zero());
	if (shapeWall->axis == 0) wallNormal = Vector3r::UnitX();
	else if (shapeWall->axis == 1)
		wallNormal = Vector3r::UnitY();
	else if (shapeWall->axis == 2)
		wallNormal = Vector3r::UnitZ();
	// Make sure that the wall normal points from the wall to the particle center
	wallNormal = (wallPosAxis - lsPosAxis > 0 ? -1 : 1) * wallNormal;

	// Get the plane equation for the wall in local coordinates
	// The plane is equation is ax + by + cz + k = 0. The normal is (a,b,c).
	Vector3r ptWallLoc = rotConjLS * (wallPos - centrEndLS);
	Vector3r nWallLoc  = rotConjLS * wallNormal;
	Real     kWallLoc  = -(ptWallLoc[0] * ptWallLoc[0] + ptWallLoc[1] * ptWallLoc[1] + ptWallLoc[2] * ptWallLoc[2]);

	// Not using the bounding ellipse sometimes causes the evaluation of a point outside of the LS grid.
	// Is there a way we can avoid this? It is unavoidable with rectangle corners penetrating the wall.
	// Using the 'unbound' variants of the normal and distance function fix this, but cost more time.

	// Box reduction based on projection on ellipsoid slice.
	Quaternionr rotGuess, rotGuessConj;
	if (useAABE && shapeLS->hasAABE) {
		// These properties should be added to the level set shape and retrieved here with something like
		Vector3r prinAxes = shapeLS->axesAABE;

		// Check if the bounding spheres overlap. If not, do not create an interaction.
		Real ctwDist = abs(wallPosAxis - lsPosAxis);
		if (!c->isReal() && !force && (ctwDist > prinAxes.maxCoeff())) {
			return false;
		} else if (ctwDist > prinAxes.maxCoeff()) {
			// Mark current interaction for removal.
			c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
			return true;
		}

		// Get a sensible guess for the normal of the interaction between the particles.
		// Normal, global coordinates, from wall to particle is wallNormal
		// Normal, local coordinates, from wall to particle is nWallLoc

		// Squared principal half axes
		Vector3r axSq = Vector3r(prinAxes[0] * prinAxes[0], prinAxes[1] * prinAxes[1], prinAxes[2] * prinAxes[2]);

		// To bound the potential overlap, find the farthest point on the ellipsoid in the direction of the guessed interaction normal
		Vector3r farPtz = -ShopLSvolume::normalToPointOnEllipsoid(nWallLoc, axSq);

		// We now proceed to make a bounding box for the overlap volume between wall and particle,
		// based on the axis-aligned bounding ellipsoid. To do this, we will rotate the ellipsoid
		// in local coordinates to align with the guessed interation normal. This allows us to
		// easily slice the ellipsoid at the farthest point in the direction of the wall normal.

		// Get the rotation that moves the local guessed interaction normal onto the local z axis (and back)
		rotGuess     = Quaternionr().setFromTwoVectors(nWallLoc, Vector3r::UnitZ());
		rotGuessConj = rotGuess.conjugate();

		// Rotate farthest points aling interaction normal onto the z axis
		Vector3r farPtzRot = rotGuess * farPtz;

		// Check if the ellipsoids can touch by using the farthest z coordinates of the ellipsoids
		if ((ctwDist + farPtzRot[2]) < 0.0) {
			// The matrix describing a general ellipsoid (not restricted to the axis-aligned case)
			Matrix3r A = Vector3r(1.0 / axSq[0], 1.0 / axSq[1], 1.0 / axSq[2]).asDiagonal();

			// Rotate the ellipsoids using the tensor relation: Anew = R*A*inv(R)
			// The coordinate systems of the ellipsoids are now aliged but with a z offset of ctcDist
			Matrix3r Arot = rotGuess * A * rotGuess.inverse();

			// We now do a plane-ellipsoid intersection at the extreme coordinates as given by the
			// farthest points on the ellipsoid in direction of the guessed interaction normal.
			// The resulting extremes of the resulting ellipses will give us the x and y limits of
			// the bounding box for the potential interaction.

			// Should funtionalise this at some point.
			// The full equation of the ellipse resulting from a sliced ellipsoid at z = h:
			// 	axx x^2 + ayy y^2 + azz h^2 + 2 axy x y + 2 axz x h + 2 ayz y h = 1
			// The h are the extreme z coordinates and locations to intersect the ellipsoids at
			Real axx = Arot(0, 0), ayy = Arot(1, 1), azz = Arot(2, 2), axy = Arot(0, 1), axz = Arot(0, 2), ayz = Arot(1, 2);

			// Location on z axis for the slicing plane, determined by extreme of the OTHER ellipsoid
			// (because that is the one that penetrates).
			// NB: For the rotated ellipsoids, particle 1 is always below particle 2!!
			Real h = -ctwDist;

			// Compute x limits of ellipse
			Real xterm1 = 0.5 / (axx - axy * axy / ayy), xterm2 = 2 * axz * h - (2 * axy * ayz * h) / ayy,
			     xterm3 = sqrt(pow(-xterm2, 2) - 4 * (-axy * axy / ayy + axx) * (-1 - ayz * ayz * h * h / ayy + azz * h * h));
			Real xlim11 = xterm1 * (xterm2 - xterm3), xlim12 = xterm1 * (xterm2 + xterm3);
			Real xlim1neg, xlim1pos;
			if (xlim11 > xlim12) {
				xlim1pos = xlim11;
				xlim1neg = xlim12;
			} else {
				xlim1pos = xlim12;
				xlim1neg = xlim11;
			}

			// Compute y limits of ellipse
			Real yterm1 = 0.5 / (ayy - axy * axy / axx), yterm2 = 2 * ayz * h - (2 * axy * axz * h) / axx,
			     yterm3 = sqrt(pow(-yterm2, 2) - 4 * (-axy * axy / axx + ayy) * (-1 - axz * axz * h * h / axx + azz * h * h));
			Real ylim11 = yterm1 * (yterm2 + yterm3), ylim12 = yterm1 * (yterm2 - yterm3);
			Real ylim1neg, ylim1pos;
			if (ylim11 > ylim12) {
				ylim1pos = ylim11;
				ylim1neg = ylim12;
			} else {
				ylim1pos = ylim12;
				ylim1neg = ylim11;
			}

			// Check if the slice we took is not actually smaller than the widest slice
			// within the potential overlap volume. Check along x axis.
			Vector3r nxGuessLoc = rotGuessConj * Vector3r::UnitX();
			Vector3r farPtx     = ShopLSvolume::normalToPointOnEllipsoid(nxGuessLoc, axSq);
			Vector3r farPtxRot  = rotGuess * farPtx;
			if (farPtxRot[2] > h) { // Point is above slice, needs to be considered
				if (farPtxRot[0] > 0) {
					xlim1pos = std::max(xlim1pos, farPtxRot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, farPtxRot[0]);
				}
				if (farPtxRot[1] > 0) {
					ylim1pos = std::max(ylim1pos, farPtxRot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, farPtxRot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			// is -farPt1xrot due to inversion symmetry.
			if (-farPtxRot[2] > h) { // Point is above slice, needs to be considered
				if (-farPtxRot[0] > 0) {
					xlim1pos = std::max(xlim1pos, -farPtxRot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, -farPtxRot[0]);
				}
				if (-farPtxRot[1] > 0) {
					ylim1pos = std::max(ylim1pos, -farPtxRot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, -farPtxRot[1]);
				}
			}
			// Now check along y axis.
			Vector3r nyGuessLoc = rotGuessConj * Vector3r::UnitY();
			Vector3r farPty     = ShopLSvolume::normalToPointOnEllipsoid(nyGuessLoc, axSq);
			Vector3r farPtyRot  = rotGuess * farPty;
			if (farPtyRot[2] > h) { // Point is above slice, needs to be considered
				if (farPtyRot[0] > 0) {
					xlim1pos = std::max(xlim1pos, farPtyRot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, farPtyRot[0]);
				}
				if (farPtyRot[1] > 0) {
					ylim1pos = std::max(ylim1pos, farPtyRot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, farPtyRot[1]);
				}
			}
			// Same for the extreme point on the other side of the ellipse
			// is -farPt1xrot due to inversion symmetry.
			if (-farPtyRot[2] > h) { // Point is above slice, needs to be considered
				if (-farPtyRot[0] > 0) {
					xlim1pos = std::max(xlim1pos, -farPtyRot[0]);
				} else {
					xlim1neg = std::min(xlim1neg, -farPtyRot[0]);
				}
				if (-farPtyRot[1] > 0) {
					ylim1pos = std::max(ylim1pos, -farPtyRot[1]);
				} else {
					ylim1neg = std::min(ylim1neg, -farPtyRot[1]);
				}
			}

			// Bounding box dimensions in local rotated frame of particle 1, intersecting the limits set by each ellipsoid
			Vector3r minBoOverlapLS = Vector3r(xlim1neg, ylim1neg, farPtzRot[2]);
			Vector3r maxBoOverlapLS = Vector3r(xlim1pos, ylim1pos, h);

			// Initialise limits in global coordinates
			Vector3r minBoOverlapOld = Vector3r(
			        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
			Vector3r maxBoOverlapOld = Vector3r(
			        -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
			// Move the old global AABB to local coordinates
			for (int i = 0; i < 2; i++) { // Go through all possible corners
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						// Select one of the corners of the AABB
						Vector3r cornerGlob = Vector3r(
						        minBoOverlap[0] + i * (sizeBoOverlap[0]),
						        minBoOverlap[1] + j * (sizeBoOverlap[1]),
						        minBoOverlap[2] + k * (sizeBoOverlap[2]));
						// Move corner to local coordinates
						Vector3r cornerLoc = rotGuess * rotConjLS * (cornerGlob - centrEndLS);
						// Find maximum and minimum values of the corners
						for (int m = 0; m < 3; m++) {
							if (cornerLoc[m] > maxBoOverlapOld[m]) { maxBoOverlapOld[m] = cornerLoc[m]; }
							if (cornerLoc[m] < minBoOverlapOld[m]) { minBoOverlapOld[m] = cornerLoc[m]; }
						}
					}
				}
			}

			// Recompute the AABB by intersecting the individual bounding boxes
			// All still local!
			for (int i = 0; i < 3; i++) {
				minBoOverlapLS[i] = math::max(minBoOverlapLS[i], minBoOverlapOld[i]);
				maxBoOverlapLS[i] = math::min(maxBoOverlapLS[i], maxBoOverlapOld[i]);
			}

			// Compute the final bounding box (BB) parameters.
			// Note: Depending on how tight the ellipsoid fits around the shape,
			// the BB can be a tight fit as well. This means that a litte bit of volume that would
			// be added because of smearing is not added. Currently, it seems that this works best.
			Vector3r cntrBoOverlapLS = (minBoOverlapLS + maxBoOverlapLS) / 2;
			cntrBoOverlap            = rotLS * (rotGuessConj * cntrBoOverlapLS) + centrEndLS;
			sizeBoOverlap            = (maxBoOverlapLS - minBoOverlapLS);
			volBoOverlap             = sizeBoOverlap.prod();
		} else {
			volBoOverlap = 0.0;
		}

		// Test if there is any overlap at all. Otherwise, do not create an interaction.
		if (!c->isReal() && !force && volBoOverlap < Mathr::EPSILON) {
			return false;
		} else if (volBoOverlap < Mathr::EPSILON) {
			// Mark current interaction for removal.
			c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
			return true;
		}
	}
	// End of overlap bounding box computations

	// Determine the initial query point (center of the overlap AABB) for the volume intergration in local axes
	// Doing this avoids that we have to do a rigid mapping for every new query point upon refining the intergration mesh
	Vector3r ptLS = rotConjLS * (cntrBoOverlap - centrEndLS);

	// Precompute a number of parameters that depend on the refinement (number of layers) of the integration mesh.
	// Start at the most coarse layer and end at the finest layer. Note that the refinement step has to be rotated as the
	// orientation of the overlap AABB is different in the respective local axes.
	std::vector<ShopLSvolume::layerData> layers(nRefineOctree); //was layers
	Real RmaxBoOverlap = 0.5 * sqrt(sizeBoOverlap[0] * sizeBoOverlap[0] + sizeBoOverlap[1] * sizeBoOverlap[1] + sizeBoOverlap[2] * sizeBoOverlap[2]);

	if (useAABE) {
		for (uint h = 0; h < nRefineOctree; h++) {
			Real     a           = 1.0 / pow(2.0, h); // Is pow() the correct function to use or is there a high-precision variant?
			Vector3r refineStep  = 0.25 * a * sizeBoOverlap;
			layers[h].cellVolume = volBoOverlap / pow(8, h); // Cell volume
			layers[h].Rmax       = a * RmaxBoOverlap;        // Maximum center-to-edge distance
			int q                = 0;                        // Looping through all the octants
			for (int i = -1; i < 2; i = i + 2) {             // Ugly way to go through all possible sign combinations
				for (int j = -1; j < 2; j = j + 2) {
					for (int k = -1; k < 2; k = k + 2) {
						Vector3r refStepPerm = Vector3r(i * refineStep[0], j * refineStep[1], k * refineStep[2]);
						// Need a different rotation because the AABB is in local coordinates of particle 1
						layers[h].refineStep1[q] = rotGuessConj * refStepPerm; // Mesh refinement step grid 1
						q++;
					}
				}
			}
		}
	} else {
		for (uint h = 0; h < nRefineOctree; h++) {
			Real     a           = 1.0 / pow(2, h); // Is pow() the correct function to use or is there a high-precision variant?
			Vector3r refineStep  = 0.25 * a * sizeBoOverlap;
			layers[h].cellVolume = volBoOverlap / pow(8, h); // Cell volume
			layers[h].Rmax       = a * RmaxBoOverlap;        // Maximum center-to-edge distance
			int q                = 0;                        // Looping through all the octants
			for (int i = -1; i < 2; i = i + 2) {             // Ugly way to go through all possible sign combinations
				for (int j = -1; j < 2; j = j + 2) {
					for (int k = -1; k < 2; k = k + 2) {
						Vector3r refStepPerm     = Vector3r(i * refineStep[0], j * refineStep[1], k * refineStep[2]) + cntrBoOverlap;
						layers[h].refineStep1[q] = rotConjLS * (refStepPerm - centrEndLS) - ptLS; // Mesh refinement step grid 1
						q++;
					}
				}
			}
		}
	}

	// Recursive volume integration
	// NB: normals and centroids are not yet normalised by the total overlap volume. Also, centroids are still in local coordinates.
	uint                            layerId = 0; // Starting layer, always zero. Needed for recursion.
	ShopLSvolume::overlapRegionData overlap
	        = ShopLSvolume::recursiveVolumeIntegrationWall(shapeLS, ptLS, nWallLoc, kWallLoc, layers, layerId, smearCoeffOctree);
	Real vn = overlap.volume;

	// Test if there is any overlap at all. Otherwise, do not create an interaction.
	if (!c->isReal() && !force && vn < Mathr::EPSILON) {
		return false;
	} else if (vn <= Mathr::EPSILON) {
		// Mark current interaction for removal.
		c->geom = ShopLSvolume::volGeomPtrForLaterRemoval(state1, state2, c);
		return true;
	}

	// Centroids and normal vector pointing from wall to particle centre
	Vector3r normal
	        = (wallNormal - (rotLS * overlap.normal1).normalized()).normalized(); // Normals need not be the same for the particles so we average them;
	// The contact point is the centroid (or center of mass) of the overlap volume
	// No averaging needed here, only one centroid was computed.
	Vector3r contactPoint = rotLS * (overlap.centroid1 / vn) + centrEndLS;

	// Compute a crude estimate of the overlap distance
	Real overlapDistance = pow(vn, 1.0 / 3.0);

	// Compute the center-to-contact distance for the torque computation.
	// NB: Only works if the normal aligns with the center-to-contact vector.  We add overlapDistance/2
	// because it is later subtracted when using e.g. ElasticContactLaw (see Law2_VolumeGeom_FrictPhys_CundallStrack::go in .cpp file).
	// This variable is now unused as the torque computation is fixed to use the proper center-to-contact vector.
	Real radius = (contactPoint - centrEndLS).norm() + 0.5 * overlapDistance; // Maybe do a single ray-trace?

	c->geom = ShopLSvolume::volGeomPtr(
	        contactPoint, // Middle of overlap volume, as usual
	        vn,           // With node contact, this would be the overlap distance, but here it is the overlap volume.
	        0,            // Because it is the wall
	        radius,
	        state1,
	        state2,
	        c,
	        normal,
	        shift2);
	return true;
}

ShopLSvolume::overlapRegionData
ShopLSvolume::recursiveVolumeIntegrationWall( // Move to ShopLS.*pp or keep here? Maybe nice to have all relevant functions in this file.
        const shared_ptr<LevelSet>&   lsShape,
        const Vector3r                pt,
        const Vector3r                nWall,
        const Real                    kWall,
        const std::vector<layerData>& layers, // cell volume, Rmax, refinement steps for grid 1 and 2
        const uint                    layerId,
        const Real                    smearCoeffOctree)
{
	// A recursive or hierarchical method for integrating the overlap volume of two intersecting level set shapes.
	// The approach is similar to the Octree data structure or adaptive mesh refinement (AMR).
	// NB: Because we allow for concave geometries, the level-set normal at the contact point need not be the correct normal
	// for computing the forces, we thus compute the normal as the volume-weighted normals of the overlap volume segments.

	// Evaluate level set values at query point
	const Real lev = lsShape->distance(pt, true);
	// Shortest distance to plane, i.e. the level set value for the infinite half space that is the wall.
	const Real levWall = -math::fabs(nWall[0] * pt[0] + nWall[1] * pt[1] + nWall[2] * pt[2] + kWall);

	// Structure for output data
	// Mabye working more with pointers would improve the code performance?
	overlapRegionData overlap;

	// Get the details of the current refinement layer
	layerData curLayer = layers[layerId];
	Real      Rmax     = curLayer.Rmax;

	// NB: lev < 0 means we are inside the particle.
	// By previous definitions, the entire AABB is already inside the wall!!
	if (lev > Rmax) {
		// END MESH REFINEMENT
		// Entire region not part of the overlap volume
		return overlap;
	} else if (lev < -Rmax) {
		// END MESH REFINEMENT
		// Entire region is part of the overlap volume, add entire cell volume to overlap volume.
		// No normal to add because there is no surface.
		Real volume       = curLayer.cellVolume; // Not sure if copying the volume in a new variable is actually slower or faster.
		overlap.volume    = volume;
		overlap.centroid1 = volume * pt;
		overlap.depth1    = volume * lev;
		overlap.depth2    = volume * levWall;
		return overlap;
	} else {
		// Region has the potential to partly belong to the overlap volume
		if (layerId < layers.size() - 1) { // I find using layers.size() ugly but I stuggled finding a more elegant end condition for the recursion.
			// REFINE THE MESH
			// Generate a new query points, dividing the region in 8 parts of equal size.
			for (int i = 0; i < 8; i++) {
				Vector3r          newPt = pt + curLayer.refineStep1[i];
				overlapRegionData overlapSub
				        = recursiveVolumeIntegrationWall(lsShape, newPt, nWall, kWall, layers, layerId + 1, smearCoeffOctree);
				overlap.volume += overlapSub.volume;
				overlap.normal1 += overlapSub.normal1;
				overlap.centroid1 += overlapSub.centroid1;
				overlap.depth1 += overlapSub.depth1;
				overlap.depth2 += overlapSub.depth2;
			}
			return overlap;
		} else {
			// END MESH REFINEMENT
			// Maximum number of refinements reached, only add volume of the cell is the centre is within both particles.
			// This assumption is reasonable because for randomly planes intersecting a rectangle, the expected volume
			// fraction of the box that is overlap is >74% if both centres as inside the box and <14% otherwise.
			Real phi = std::max(lev, levWall);
			Real phiRef(0.);
			if (smearCoeffOctree != 0) phiRef = Rmax / smearCoeffOctree;
			if (abs(phi) < phiRef) {
				// Compute the smeared Heaviside step function
				Real volume       = 0.5 * (1.0 - phi / phiRef + sin(-M_PI * phi / phiRef) / M_PI) * curLayer.cellVolume;
				overlap.volume    = volume;
				overlap.centroid1 = volume * pt;
				overlap.depth1    = volume * lev;
				overlap.depth2    = volume * levWall;
				Real surface      = 0.5 / phiRef * (1 + cos(M_PI * phi / phiRef)) * pow(curLayer.cellVolume, 2.0 / 3.0);
				overlap.area      = surface;
				overlap.normal1   = surface * (lsShape->normal(pt, true));
				return overlap;
			} else if (phi < 0) {
				// No normal to add because there is no surface.
				Real volume       = curLayer.cellVolume;
				overlap.volume    = volume;
				overlap.centroid1 = volume * pt;
				overlap.depth1    = volume * lev;
				overlap.depth2    = volume * levWall;
				return overlap;
			} else {
				return overlap;
			}
		}
	} // End mesh refinement check
}

} // namespace yade
#endif //YADE_LS_DEM
