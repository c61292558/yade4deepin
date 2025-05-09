
#ifdef FLOW_ENGINE

#if CGAL_VERSION_NR < CGAL_VERSION_NUMBER(4, 11, 0)
#include "CGAL/constructions/constructions_on_weighted_points_cartesian_3.h"
#endif
#include "vector"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <new>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>


namespace yade { // Cannot have #include directive inside.

namespace CGT {


	// 	template<class Tesselation> const Real Network<Tesselation>::FAR = 50000;
	template <class Tesselation> const Real Network<Tesselation>::ONE_THIRD = 1.0 / 3.0;
	template <class Tesselation> const int  Network<Tesselation>::facetVertices[4][3] = { { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };
	template <class Tesselation> const int  Network<Tesselation>::permut3[3][3] = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } };
	template <class Tesselation> const int  Network<Tesselation>::permut4[4][4] = { { 0, 1, 2, 3 }, { 1, 2, 3, 0 }, { 2, 3, 0, 1 }, { 3, 0, 1, 2 } };

	template <class Tesselation> Network<Tesselation>::~Network() { }

	template <class Tesselation> Network<Tesselation>::Network()
	{
		FAR = 50000;
		facetF1 = facetF2 = facetRe1 = facetRe2 = facetRe3 = 0;
		// 	F1=F2=Re1=Re2=0;
	}

	template <class Tesselation> int Network<Tesselation>::detectFacetFictiousVertices(CellHandle& cell, int& j)
	{
		facetNFictious = 0;
		int nRealVtx = 0;
		for (int kk = 0; kk < 3; kk++) {
			if (cell->vertex(facetVertices[j][kk])->info().isFictious) {
				if (facetNFictious == 0) facetF1 = kk;
				else
					facetF2 = kk;
				facetNFictious += 1;
			} else {
				if (nRealVtx == 0) facetRe1 = kk;
				else if (nRealVtx == 1)
					facetRe2 = kk;
				else if (nRealVtx == 2)
					facetRe3 = kk;
				nRealVtx += 1;
			}
		}
		return facetNFictious;
	}

	template <class Tesselation> Real Network<Tesselation>::volumePoreVoronoiFraction(CellHandle& cell, int& j, bool reuseFacetData)
	{
		Point& p1 = cell->info();
		Point& p2 = cell->neighbor(j)->info();
		if (!reuseFacetData) facetNFictious = detectFacetFictiousVertices(cell, j);
		//cout << "getting volume pore voronoi fraction for alpha ? "<< cell->info().isAlpha << " fictiousN " << facetNFictious << endl;
		Sphere       v[3];
		VertexHandle W[3];
		for (int kk = 0; kk < 3; kk++) {
			W[kk] = cell->vertex(facetVertices[j][kk]);
			v[kk] = cell->vertex(facetVertices[j][kk])->point();
		}
		switch (facetNFictious) {
			case (0): {
				VertexHandle& SV1 = W[0];
				VertexHandle& SV2 = W[1];
				VertexHandle& SV3 = W[2];

				cell->info().facetSurfaces[j]
				        = 0.5 * CGAL::cross_product(SV1->point().point() - SV3->point().point(), SV2->point().point() - SV3->point().point());
				if (cell->info().facetSurfaces[j][0] == 0 && cell->info().facetSurfaces[j][1] == 0 && cell->info().facetSurfaces[j][2] == 0)
					cerr << "NULL FACET SURF" << endl;
				if (cell->info().facetSurfaces[j] * (p2 - p1) > 0) cell->info().facetSurfaces[j] = -1.0 * cell->info().facetSurfaces[j];
				Real Vtot = abs(ONE_THIRD * cell->info().facetSurfaces[j] * (p1 - p2));
				Vtotalissimo += Vtot;

				Real Vsolid1 = 0, Vsolid2 = 0;
				for (int i = 0; i < 3; i++) {
					Vsolid1 += sphericalTriangleVolume(v[permut3[i][0]], v[permut3[i][1]].point(), p1, p2);
					Vsolid2 += sphericalTriangleVolume(v[permut3[i][0]], v[permut3[i][2]].point(), p1, p2);
				}

				VSolidTot += Vsolid1 + Vsolid2;
				vPoral += Vtot - (Vsolid1 + Vsolid2);

				bool border = false;
				for (int i = 0; i < 4; i++) {
					if (cell->neighbor(i)->info().fictious() != 0) border = true;
				}
				if (!border) {
					vPoralPorosity += Vtot - (Vsolid1 + Vsolid2);
					vTotalPorosity += Vtot;
				}

				/**Vpore**/ return Vtot - (Vsolid1 + Vsolid2);
			}; break;
			case (1): {
				return volumeSingleFictiousPore(
				        cell->vertex(facetVertices[j][facetF1]),
				        cell->vertex(facetVertices[j][facetRe1]),
				        cell->vertex(facetVertices[j][facetRe2]),
				        p1,
				        p2,
				        cell->info().facetSurfaces[j]);
			}; break;
			case (2): {
				return volumeDoubleFictiousPore(
				        cell->vertex(facetVertices[j][facetF1]),
				        cell->vertex(facetVertices[j][facetF2]),
				        cell->vertex(facetVertices[j][facetRe1]),
				        p1,
				        p2,
				        cell->info().facetSurfaces[j]);
			}; break;
			default: return 0;
		}
	}

	template <class Tesselation> Real Network<Tesselation>::volumeSolidPore(const CellHandle& cell)
	{
		Real Vsolid = 0;
		for (int i = 0; i < 4; i++) {
			if (!cell->vertex(permut4[i][0])->info().isFictious)
				Vsolid += sphericalTriangleVolume(
				        cell->vertex(permut4[i][0])->point(),
				        cell->vertex(permut4[i][1])->point().point(),
				        cell->vertex(permut4[i][2])->point().point(),
				        cell->vertex(permut4[i][3])->point().point());
		}
		return Vsolid;
	}

	template <class Tesselation>
	Real Network<Tesselation>::volumeSingleFictiousPore(
	        const VertexHandle& SV1, const VertexHandle& SV2, const VertexHandle& SV3, const Point& PV1, const Point& PV2, CVector& facetSurface)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		Real      A[3], B[3];
		Boundary& bi1 = boundary(SV1->info().id());

		for (int m = 0; m < 3; m++) {
			A[m] = (SV2->point())[m];
		}
		for (int m = 0; m < 3; m++) {
			B[m] = (SV3->point())[m];
		}

		A[bi1.coordinate] = bi1.p[bi1.coordinate];
		B[bi1.coordinate] = bi1.p[bi1.coordinate];

		Point AA(A[0], A[1], A[2]);
		Point BB(B[0], B[1], B[2]);
		facetSurface = surfaceSingleFictiousFacet(SV1, SV2, SV3);
		if (facetSurface * (PV2 - PV1) > 0) facetSurface = -1.0 * facetSurface;
		Real Vtot = ONE_THIRD * abs(facetSurface * (PV1 - PV2));
		Vtotalissimo += Vtot;

		Sphere  A1(AA, 0);
		Sphere  B1(BB, 0);
		Sphere& SW2 = SV2->point();
		Sphere& SW3 = SV3->point();

		Real Vsolid1 = sphericalTriangleVolume(SW2, AA, PV1, PV2) + sphericalTriangleVolume(SW2, SW3.point(), PV1, PV2);
		Real Vsolid2 = sphericalTriangleVolume(SW3, BB, PV1, PV2) + sphericalTriangleVolume(SW3, SW2.point(), PV1, PV2);

		VSolidTot += Vsolid1 + Vsolid2;
		vPoral += Vtot - (Vsolid1 + Vsolid2);

		return (Vtot - (Vsolid1 + Vsolid2));
	}

	template <class Tesselation>
	Real Network<Tesselation>::volumeDoubleFictiousPore(
	        const VertexHandle& SV1, const VertexHandle& SV2, const VertexHandle& SV3, const Point& PV1, const Point& PV2, CVector& facetSurface)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		Real      A[3], B[3];
		Boundary& bi1 = boundary(SV1->info().id());
		Boundary& bi2 = boundary(SV2->info().id());
		for (int m = 0; m < 3; m++) {
			A[m] = B[m] = SV3->point()[m];
		}

		A[bi1.coordinate] = bi1.p[bi1.coordinate];
		B[bi2.coordinate] = bi2.p[bi2.coordinate];
		Point AA(A[0], A[1], A[2]);
		Point BB(B[0], B[1], B[2]);

		facetSurface = CGAL::cross_product(SV3->point().point() - AA, SV3->point().point() - BB);
		if (facetSurface * (PV2 - PV1) > 0) facetSurface = -1.0 * facetSurface;
		Real Vtot = abs(facetSurface * (PV1 - PV2)) * ONE_THIRD;
		Vtotalissimo += Vtot;

		Real Vsolid1 = sphericalTriangleVolume(SV3->point(), AA, PV1, PV2);
		Real Vsolid2 = sphericalTriangleVolume(SV3->point(), BB, PV1, PV2);

		vPoral += (Vtot - Vsolid1 - Vsolid2);
		VSolidTot += Vsolid1 + Vsolid2;

		return (Vtot - Vsolid1 - Vsolid2);
	}

	template <class Tesselation> Real Network<Tesselation>::sphericalTriangleVolume(const Sphere& ST1, const Point& PT1, const Point& PT2, const Point& PT3)
	{
		Real rayon = sqrt(ST1.weight());
		if (rayon == 0.0) return 0.0;
		return ((ONE_THIRD * rayon) * (fastSphericalTriangleArea(ST1, PT1, PT2, PT3)));
	}

	template <class Tesselation>
	Real Network<Tesselation>::fastSphericalTriangleArea(const Sphere& STA1, const Point& STA2, const Point& STA3, const Point& PTA1)
	{
		using namespace CGAL;
		Real rayon2 = STA1.weight();
		if (rayon2 == 0.0) return 0.0;
		return rayon2 * fastSolidAngle(STA1.point(), STA2, STA3, PTA1);
	}

	template <class Tesselation> Real Network<Tesselation>::sphericalTriangleArea(Sphere STA1, Sphere STA2, Sphere STA3, Point PTA1)
	{
		Real rayon = STA1.weight();
		if (rayon == 0.0) return 0.0;

		CVector v12 = STA2.point() - STA1.point();
		CVector v13 = STA3.point() - STA1.point();
		CVector v14 = PTA1 - STA1.point();

		Real norme12 = (v12.squared_length());
		Real norme13 = (v13.squared_length());
		Real norme14 = (v14.squared_length());

		Real cosA = v12 * v13 / sqrt((norme13 * norme12));
		Real cosB = v12 * v14 / sqrt((norme14 * norme12));
		Real cosC = v14 * v13 / sqrt((norme13 * norme14));

		Real A = acos(cosA);
		Real B = acos(cosB);
		Real C = acos(cosC);
		if (A == 0 || B == 0 || C == 0) return 0;

		Real a = acos((cosA - cosB * cosC) / (sin(B) * sin(C)));
		Real b = acos((cosB - cosC * cosA) / (sin(C) * sin(A)));
		Real c = acos((cosC - cosA * cosB) / (sin(A) * sin(B)));

		Real aire_triangle_spherique = rayon * (a + b + c - M_PI);

		return aire_triangle_spherique;
	}

	template <class Tesselation> Real Network<Tesselation>::fastSolidAngle(const Point& STA1, const Point& PTA1, const Point& PTA2, const Point& PTA3)
	{
		//! This function needs to be fast because it is used heavily. Avoid using vector operations which require constructing vectors (~50% of cpu time in the non-fast version), and compute angle using the 3x faster formula of Oosterom and StrackeeVan Oosterom, A; Strackee, J (1983). "The Solid Angle of a Plane Triangle". IEEE Trans. Biom. Eng. BME-30 (2): 125-126. (or check http://en.wikipedia.org/wiki/Solid_angle)
		using namespace CGAL;
		Real M[3][3];
		M[0][0] = PTA1.x() - STA1.x();
		M[0][1] = PTA2.x() - STA1.x();
		M[0][2] = PTA3.x() - STA1.x();
		M[1][0] = PTA1.y() - STA1.y();
		M[1][1] = PTA2.y() - STA1.y();
		M[1][2] = PTA3.y() - STA1.y();
		M[2][0] = PTA1.z() - STA1.z();
		M[2][1] = PTA2.z() - STA1.z();
		M[2][2] = PTA3.z() - STA1.z();

		Real detM = M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) + M[1][0] * (M[2][1] * M[0][2] - M[0][1] * M[2][2])
		        + M[2][0] * (M[0][1] * M[1][2] - M[1][1] * M[0][2]);

		Real pv12N2 = pow(M[0][0], 2) + pow(M[1][0], 2) + pow(M[2][0], 2);
		Real pv13N2 = pow(M[0][1], 2) + pow(M[1][1], 2) + pow(M[2][1], 2);
		Real pv14N2 = pow(M[0][2], 2) + pow(M[1][2], 2) + pow(M[2][2], 2);

		Real pv12N = sqrt(pv12N2);
		Real pv13N = sqrt(pv13N2);
		Real pv14N = sqrt(pv14N2);

		Real cp12 = (M[0][0] * M[0][1] + M[1][0] * M[1][1] + M[2][0] * M[2][1]);
		Real cp13 = (M[0][0] * M[0][2] + M[1][0] * M[1][2] + M[2][0] * M[2][2]);
		Real cp23 = (M[0][1] * M[0][2] + M[1][1] * M[1][2] + M[2][1] * M[2][2]);

		Real ratio = detM / (pv12N * pv13N * pv14N + cp12 * pv14N + cp13 * pv13N + cp23 * pv12N);
		return math::abs(2 * atan(ratio));
	}

	template <class Tesselation> Real Network<Tesselation>::surfaceSolidThroat(CellHandle cell, int j, bool slipBoundary, bool reuseFacetData)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		if (!reuseFacetData) facetNFictious = detectFacetFictiousVertices(cell, j);
		Point& p1 = cell->info();
		Point& p2 = cell->neighbor(j)->info();

		Real Ssolid = 0;
		Real Ssolid1 = 0, Ssolid1n = 0, Ssolid2 = 0, Ssolid2n = 0, Ssolid3 = 0, Ssolid3n = 0;

		Sphere       v[3];
		VertexHandle W[3];

		for (int kk = 0; kk < 3; kk++) {
			W[kk] = cell->vertex(facetVertices[j][kk]);
			v[kk] = cell->vertex(facetVertices[j][kk])->point();
		}

		switch (facetNFictious) {
			case (0): {
				VertexHandle& SV1 = W[0];
				VertexHandle& SV2 = W[1];
				VertexHandle& SV3 = W[2];

				Ssolid1 = fastSphericalTriangleArea(SV1->point(), SV2->point().point(), p1, p2);
				Ssolid1n = fastSphericalTriangleArea(SV1->point(), SV3->point().point(), p1, p2);
				cell->info().solidSurfaces[j][0] = Ssolid1 + Ssolid1n;
				Ssolid2 = fastSphericalTriangleArea(SV2->point(), SV1->point().point(), p1, p2);
				Ssolid2n = fastSphericalTriangleArea(SV2->point(), SV3->point().point(), p1, p2);
				cell->info().solidSurfaces[j][1] = Ssolid2 + Ssolid2n;
				Ssolid3 = fastSphericalTriangleArea(SV3->point(), SV2->point().point(), p1, p2);
				Ssolid3n = fastSphericalTriangleArea(SV3->point(), SV1->point().point(), p1, p2);
				cell->info().solidSurfaces[j][2] = Ssolid3 + Ssolid3n;

			}; break;
			case (1): {
				VertexHandle SV1 = cell->vertex(facetVertices[j][facetF1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetRe1]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetRe2]);

				Boundary& bi1 = boundary(SV1->info().id());
				Ssolid1 = 0;
				if (bi1.flowCondition && !slipBoundary) {
					Ssolid1 = abs(0.5 * CGAL::cross_product(p1 - p2, SV2->point().point() - SV3->point().point())[bi1.coordinate]);
					cell->info().solidSurfaces[j][facetF1] = Ssolid1;
				}
				Ssolid2 = fastSphericalTriangleArea(SV2->point(), SV1->point().point(), p1, p2);
				Ssolid2n = fastSphericalTriangleArea(SV2->point(), SV3->point().point(), p1, p2);
				cell->info().solidSurfaces[j][facetRe1] = Ssolid2 + Ssolid2n;
				Ssolid3 = fastSphericalTriangleArea(SV3->point(), SV2->point().point(), p1, p2);
				Ssolid3n = fastSphericalTriangleArea(SV3->point(), SV1->point().point(), p1, p2);
				cell->info().solidSurfaces[j][facetRe2] = Ssolid3 + Ssolid3n;
			}; break;
			case (2): {
				Real A[3], B[3], C[3];

				VertexHandle SV1 = cell->vertex(facetVertices[j][facetF1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetF2]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetRe1]);

				Boundary& bi1 = boundary(SV1->info().id());
				Boundary& bi2 = boundary(SV2->info().id());
				for (int m = 0; m < 3; m++) {
					A[m] = B[m] = C[m] = (SV3->point())[m];
				}
				A[bi1.coordinate] = bi1.p[bi1.coordinate];
				B[bi2.coordinate] = bi2.p[bi2.coordinate];
				C[bi1.coordinate] = bi1.p[bi1.coordinate];
				C[bi2.coordinate] = bi2.p[bi2.coordinate];
				Point AA(A[0], A[1], A[2]);
				Point BB(B[0], B[1], B[2]);
				Point CC(C[0], C[1], C[2]);

				Sphere A1(AA, 0);
				Sphere B1(BB, 0);
				Sphere C1(CC, 0);
				//FIXME : we are computing triangle area twice here, because its computed in volume_double_fictious already -> optimize
				Ssolid1 = fastSphericalTriangleArea(SV3->point(), AA, p1, p2);
				Ssolid1n = fastSphericalTriangleArea(SV3->point(), BB, p1, p2);
				cell->info().solidSurfaces[j][facetRe1] = Ssolid1 + Ssolid1n;
				//area vector of triangle (p1,sphere,p2)
				CVector p1p2v1Surface = 0.5 * CGAL::cross_product(p1 - p2, SV3->point().point() - p2);
				if (bi1.flowCondition && !slipBoundary) {
					//projection on boundary 1
					Ssolid2 = abs(p1p2v1Surface[bi1.coordinate]);
					cell->info().solidSurfaces[j][facetF1] = Ssolid2;
				} else
					cell->info().solidSurfaces[j][facetF1] = 0;

				if (bi2.flowCondition && !slipBoundary) {
					//projection on boundary 2
					Ssolid3 = abs(p1p2v1Surface[bi2.coordinate]);
					cell->info().solidSurfaces[j][facetF2] = Ssolid3;
				} else
					cell->info().solidSurfaces[j][facetF2] = 0;
			}; break;
			default: throw std::runtime_error(__FILE__ " : switch default case error.");
		}

		Ssolid = Ssolid1 + Ssolid1n + Ssolid2 + Ssolid2n + Ssolid3 + Ssolid3n;

		if (Ssolid) cell->info().solidSurfaces[j][3] = 1.0 / Ssolid;
		else
			cell->info().solidSurfaces[j][3] = 0;
		sSolidTot += Ssolid;

		return Ssolid;
	}

	template <class Tesselation> Real Network<Tesselation>::surfaceSolidThroatInPore(CellHandle cell, int j, bool slipBoundary, bool reuseFacetData)
	{
		using math::abs; // inside function it does not leak.
		if (!reuseFacetData) facetNFictious = detectFacetFictiousVertices(cell, j);
		Point&       p1 = cell->info();
		Point&       p2 = cell->neighbor(j)->info();
		Real         Ssolid1 = 0, Ssolid2 = 0, Ssolid3 = 0;
		Sphere       v[3];
		VertexHandle W[3];

		for (int kk = 0; kk < 3; kk++) {
			W[kk] = cell->vertex(facetVertices[j][kk]);
			v[kk] = cell->vertex(facetVertices[j][kk])->point();
		}

		switch (facetNFictious) {
			case (0): {
				VertexHandle& SV1 = W[0];
				VertexHandle& SV2 = W[1];
				VertexHandle& SV3 = W[2];

				Ssolid1 = fastSphericalTriangleArea(SV1->point(), SV2->point().point(), p1, SV3->point().point());
				Ssolid2 = fastSphericalTriangleArea(SV2->point(), SV1->point().point(), p1, SV3->point().point());
				Ssolid3 = fastSphericalTriangleArea(SV3->point(), SV2->point().point(), p1, SV1->point().point());
			}; break;
			case (1): {
				VertexHandle SV1 = cell->vertex(facetVertices[j][facetF1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetRe1]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetRe2]);

				Boundary& bi1 = boundary(SV1->info().id());
				Ssolid1 = 0;
				if (bi1.flowCondition && !slipBoundary)
					Ssolid1 = abs(
					        0.5
					        * CGAL::cross_product(p1 - SV2->point().point(), SV2->point().point() - SV3->point().point())[bi1.coordinate]);
				Ssolid2 = fastSphericalTriangleArea(SV2->point(), SV3->point().point(), p1, SV2->point().point() + bi1.normal);
				Ssolid3 = fastSphericalTriangleArea(SV3->point(), SV2->point().point(), p1, SV3->point().point() + bi1.normal);
			}; break;
			case (2): {
				Real         A[3], B[3], C[3];
				VertexHandle SV1 = cell->vertex(facetVertices[j][facetF1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetF2]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetRe1]);
				Boundary&    bi1 = boundary(SV1->info().id());
				Boundary&    bi2 = boundary(SV2->info().id());

				for (int m = 0; m < 3; m++) {
					A[m] = B[m] = C[m] = (SV3->point())[m];
				}
				A[bi1.coordinate] = bi1.p[bi1.coordinate];
				B[bi2.coordinate] = bi2.p[bi2.coordinate];
				C[bi1.coordinate] = bi1.p[bi1.coordinate];
				C[bi2.coordinate] = bi2.p[bi2.coordinate];
				Point AA(A[0], A[1], A[2]);
				Point BB(B[0], B[1], B[2]);
				Point CC(C[0], C[1], C[2]);

				Sphere A1(AA, 0);
				Sphere B1(BB, 0);
				Sphere C1(CC, 0);
				//FIXME : we are computing triangle area twice here, because its computed in volume_double_fictious already -> optimize
				Ssolid1 = 0.5 * (fastSphericalTriangleArea(SV3->point(), AA, p1, p2) + fastSphericalTriangleArea(SV3->point(), BB, p1, p2));
				CVector p1p2v1Surface = 0.5 * CGAL::cross_product(p1 - p2, SV3->point().point() - p2);
				if (bi1.flowCondition && !slipBoundary) Ssolid2 = 0.5 * abs(p1p2v1Surface[bi1.coordinate]);
				if (bi2.flowCondition && !slipBoundary) Ssolid3 = 0.5 * abs(p1p2v1Surface[bi2.coordinate]);
			}; break;
			default: throw std::runtime_error(__FILE__ " : switch default case error.");
		}
		return Ssolid1 + Ssolid2 + Ssolid3;
	}

	template <class Tesselation> CVector Network<Tesselation>::surfaceDoubleFictiousFacet(VertexHandle fSV1, VertexHandle fSV2, VertexHandle SV3)
	{
		//This function is correct only with axis-aligned boundaries FIXME: alpha is not axis aligned, is there a problem?
		const Boundary& bi1 = boundary(fSV1->info().id());
		const Boundary& bi2 = boundary(fSV2->info().id());

		Real area = (bi1.p[bi1.coordinate] - SV3->point()[bi1.coordinate]) * (bi2.p[bi2.coordinate] - SV3->point()[bi2.coordinate]);
		Real surf[3] = { 1, 1, 1 };
		surf[bi1.coordinate] = 0;
		surf[bi2.coordinate] = 0;
		return area * CVector(surf[0], surf[1], surf[2]);
	}

	template <class Tesselation> CVector Network<Tesselation>::surfaceSingleFictiousFacet(VertexHandle fSV1, VertexHandle SV2, VertexHandle SV3)
	{
		//This function is correct only with axis-aligned boundaries
		const Boundary& bi1 = boundary(fSV1->info().id());
		//  const Boundary &bi2 = boundary ( fSV2->info().id() );
		CVector mean_height = (bi1.p[bi1.coordinate] - 0.5 * (SV3->point()[bi1.coordinate] + SV2->point()[bi1.coordinate])) * bi1.normal;

		return CGAL::cross_product(mean_height, SV3->point().point() - SV2->point().point());
	}

	template <class Tesselation> Real Network<Tesselation>::surfaceSolidDoubleFictiousFacet(VertexHandle SV1, VertexHandle SV2, VertexHandle SV3)
	{
		Real A[3], B[3];

		for (int m = 0; m < 3; m++) {
			A[m] = B[m] = (SV3->point())[m];
		}

		const Boundary& bi1 = boundary(SV1->info().id());
		const Boundary& bi2 = boundary(SV2->info().id());

		A[bi1.coordinate] = bi1.p[bi1.coordinate];
		B[bi2.coordinate] = bi2.p[bi2.coordinate];

		Real board1 = 0, board2 = 0, total_surface = 0;
		for (int p = 0; p < 3; p++) {
			board1 += pow((SV3->point()[p] - A[p]), 2);
			board2 += pow((SV3->point()[p] - B[p]), 2);
		}
		total_surface = sqrt(board1 * board2);

		Real solid_surface = surfaceSolidFacet(SV3->point(), SV2->point(), SV1->point());

		return total_surface - solid_surface;
	}

	template <class Tesselation> Real Network<Tesselation>::surfaceSolidFacet(Sphere ST1, Sphere ST2, Sphere ST3)
	{
		Real Area;

		Real squared_radius = ST1.weight();

		CVector v12 = ST2.point() - ST1.point();
		CVector v13 = ST3.point() - ST1.point();

		Real norme12 = v12.squared_length();
		Real norme13 = v13.squared_length();

		Real cosA = v12 * v13 / (sqrt(norme13 * norme12));
		Real A = acos(cosA);

		Area = (A / (2 * M_PI)) * (M_PI * squared_radius);
		return Area;
	}

	template <class Tesselation> void Network<Tesselation>::addBoundingPlanes()
	{
		Tesselation& Tes = T[currentTes];
		//FIXME: Id's order in boundsIds is done according to the enumerotation of boundaries from TXStressController.hpp, line 31. DON'T CHANGE IT!
		yMinId = Tes.Max_id() + 2;
		boundsIds[0] = &yMinId;
		yMaxId = Tes.Max_id() + 3;
		boundsIds[1] = &yMaxId;
		xMinId = Tes.Max_id() + 0;
		boundsIds[2] = &xMinId;
		xMaxId = Tes.Max_id() + 1;
		boundsIds[3] = &xMaxId;
		zMinId = Tes.Max_id() + 5;
		boundsIds[4] = &zMaxId;
		zMaxId = Tes.Max_id() + 6;
		boundsIds[5] = &zMinId;

		cornerMin = Point(xMin, yMin, zMin);
		cornerMax = Point(xMax, yMax, zMax);

		idOffset
		        = Tes.Max_id() + 1; //so that boundaries[vertex->id - offset] gives the ordered boundaries (also see function Boundary& boundary(int b))

		addBoundingPlane(CVector(0, 1, 0), yMinId);
		addBoundingPlane(CVector(0, -1, 0), yMaxId);
		addBoundingPlane(CVector(-1, 0, 0), xMaxId);
		addBoundingPlane(CVector(1, 0, 0), xMinId);
		addBoundingPlane(CVector(0, 0, 1), zMinId);
		addBoundingPlane(CVector(0, 0, -1), zMaxId);

		// 	addBoundingPlanes(true);
	}

	template <class Tesselation> void Network<Tesselation>::addBoundingPlane(CVector Normal, int id_wall)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		                 // 	  Tesselation& Tes = T[currentTes];
		//FIXME: pre-condition: the normal is axis-aligned
		int Coordinate
		        = int(math::round(math::abs(Normal[0]))) * 0 + int(math::round(math::abs(Normal[1]))) * 1 + int(math::round(math::abs(Normal[2]))) * 2;

		Real pivot = Normal[Coordinate] < 0 ? cornerMax.x() * abs(Normal[0]) + cornerMax.y() * abs(Normal[1]) + cornerMax.z() * abs(Normal[2])
		                                    : cornerMin.x() * abs(Normal[0]) + cornerMin.y() * abs(Normal[1]) + cornerMin.z() * abs(Normal[2]);

		Real center[3] = { 0.5 * (cornerMin.x() + cornerMax.x()) * (1 - abs(Normal[0])) + pivot * abs(Normal[0]),
			           0.5 * (cornerMax.y() + cornerMin.y()) * (1 - abs(Normal[1])) + pivot * abs(Normal[1]),
			           0.5 * (cornerMax.z() + cornerMin.z()) * (1 - abs(Normal[2])) + pivot * abs(Normal[2]) };

		addBoundingPlane(center, 0, Normal, id_wall);
	}

	template <class Tesselation> void Network<Tesselation>::addBoundingPlane(Real center[3], Real thickness, CVector Normal, int id_wall)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		Tesselation& Tes = T[currentTes];

		int Coordinate
		        = int(math::round(math::abs(Normal[0]))) * 0 + int(math::round(math::abs(Normal[1]))) * 1 + int(math::round(math::abs(Normal[2]))) * 2;

		Tes.insert(
		        (center[0] + Normal[0] * thickness / 2) * (1 - abs(Normal[0]))
		                + (center[0] + Normal[0] * thickness / 2 - Normal[0] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[0]),
		        (center[1] + Normal[1] * thickness / 2) * (1 - abs(Normal[1]))
		                + (center[1] + Normal[1] * thickness / 2 - Normal[1] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[1]),
		        (center[2] + Normal[2] * thickness / 2) * (1 - abs(Normal[2]))
		                + (center[2] + Normal[2] * thickness / 2 - Normal[2] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[2]),
		        FAR * (cornerMax.y() - cornerMin.y()),
		        id_wall,
		        true);

		Point P(center[0], center[1], center[2]);
		boundaries[id_wall - idOffset].p = P;
		boundaries[id_wall - idOffset].normal = Normal;
		boundaries[id_wall - idOffset].coordinate = Coordinate;

		boundaries[id_wall - idOffset].flowCondition = 1;
		boundaries[id_wall - idOffset].value = 0;

		if (debugOut)
			cout << "A boundary -center/thick- has been created. ID = " << id_wall << " position = "
			     << (center[0] + Normal[0] * thickness / 2) * (1 - abs(Normal[0]))
			                + (center[0] + Normal[0] * thickness / 2 - Normal[0] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[0])
			     << " , "
			     << (center[1] + Normal[1] * thickness / 2) * (1 - abs(Normal[1]))
			                + (center[1] + Normal[1] * thickness / 2 - Normal[1] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[1])
			     << " , "
			     << (center[2] + Normal[2] * thickness / 2) * (1 - abs(Normal[2]))
			                + (center[2] + Normal[2] * thickness / 2 - Normal[2] * FAR * (cornerMax.y() - cornerMin.y())) * abs(Normal[2])
			     << ". Radius = " << FAR * (cornerMax.y() - cornerMin.y()) << endl;
	}

	template <class Tesselation> void Network<Tesselation>::defineFictiousCells()
	{
		RTriangulation&     Tri = T[currentTes].Triangulation();
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			cell->info().fictious() = 0;
		}
		for (int bound = 0; bound < 6; bound++) {
			int& id = *boundsIds[bound];
			if (id < 0) continue;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunused-result"
			T[currentTes].vertexHandles[id];
#pragma GCC diagnostic pop
			VectorCell tmpCells;
			tmpCells.resize(10000); // limiting if we have >10000 fictitious cells?
			VCellIterator cells_it = tmpCells.begin();
			VCellIterator cells_end = Tri.incident_cells(T[currentTes].vertexHandles[id], cells_it); //
			for (VCellIterator it = tmpCells.begin(); it != cells_end; it++) {
				CellHandle& cell = *it;
				(cell->info().fictious()) += 1;
				cell->info().isFictious = true;
			}
		}
		if (debugOut) cout << "Fictious cell defined" << endl;
	}


	// Real Network::sphericalTriangleArea ( Sphere STA1, Sphere STA2, Sphere STA3, Point PTA1 )
	// {
	//  Real rayon = STA1.weight();
	//  if ( rayon == 0.0 ) return 0.0;
	//
	//  CVector v12 = STA2.point() - STA1.point();
	//  CVector v13 = STA3.point() - STA1.point();
	//  CVector v14 = PTA1 - STA1.point();
	//
	//  Real norme12 = ( v12.squared_length() );
	//  Real norme13 = ( v13.squared_length() );
	//  Real norme14 = ( v14.squared_length() );
	//
	//  Real cosA = v12*v13 / sqrt ( ( norme13 * norme12 ) );
	//  Real cosB = v12*v14 / sqrt ( ( norme14 * norme12 ) );
	//  Real cosC = v14*v13 / sqrt ( ( norme13 * norme14 ) );
	//
	//  Real A = acos ( cosA );
	//  Real B = acos ( cosB );
	//  Real C = acos ( cosC );
	//  if ( A==0 || B==0 || C==0 ) return 0;
	//
	//  Real a = acos ( ( cosA - cosB * cosC ) / ( sin ( B ) * sin ( C ) ) );
	//  Real b = acos ( ( cosB - cosC * cosA ) / ( sin ( C ) * sin ( A ) ) );
	//  Real c = acos ( ( cosC - cosA * cosB ) / ( sin ( A ) * sin ( B ) ) );
	//
	//  Real aire_triangle_spherique = rayon * ( a + b + c - M_PI );
	//
	//  return  aire_triangle_spherique;
	// }

	template <class Tesselation> void Network<Tesselation>::lineSolidPore(CellHandle cell, int j)
	{
		using math::abs; // when used inside function it does not leak - it is safe.
		facetNFictious = detectFacetFictiousVertices(cell, j);
		Real         lSolid = 0; //total of solidLine[j][0], solidLine[j][1], solidLine[j][2].
		Sphere       v[3];
		VertexHandle W[3];

		for (int kk = 0; kk < 3; kk++) {
			W[kk] = cell->vertex(facetVertices[j][kk]);
			v[kk] = cell->vertex(facetVertices[j][kk])->point();
		}

		switch (facetNFictious) {
			case (0): {
				VertexHandle& SV1 = W[0];
				VertexHandle& SV2 = W[1];
				VertexHandle& SV3 = W[2];

				cell->info().solidLine[j][0] = lineSolidFacet(SV1->point(), SV2->point(), SV3->point());
				cell->info().solidLine[j][1] = lineSolidFacet(SV2->point(), SV3->point(), SV1->point());
				cell->info().solidLine[j][2] = lineSolidFacet(SV3->point(), SV1->point(), SV2->point());
			}; break;
			case (1): {
				VertexHandle SV1 = cell->vertex(facetVertices[j][facetRe1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetRe2]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetF1]);

				cell->info().solidLine[j][facetRe1] = lineSolidFacet(SV1->point(), SV2->point(), SV3->point());
				cell->info().solidLine[j][facetRe2] = lineSolidFacet(SV2->point(), SV1->point(), SV3->point());

				Boundary& bi = boundary(SV3->info().id());
				Real      A[3], B[3];
				for (int m = 0; m < 3; m++) {
					A[m] = SV1->point()[m];
					B[m] = SV2->point()[m];
				}
				A[bi.coordinate] = bi.p[bi.coordinate];
				B[bi.coordinate] = bi.p[bi.coordinate];
				Point   AA(A[0], A[1], A[2]);
				Point   BB(B[0], B[1], B[2]);
				CVector AB = AA - BB;
				cell->info().solidLine[j][facetF1] = sqrt(AB.squared_length());
			}; break;
			case (2): {
				VertexHandle SV1 = cell->vertex(facetVertices[j][facetF1]);
				VertexHandle SV2 = cell->vertex(facetVertices[j][facetF2]);
				VertexHandle SV3 = cell->vertex(facetVertices[j][facetRe1]);

				cell->info().solidLine[j][facetRe1] = 0.5 * M_PI * sqrt(SV3->point().weight());

				Boundary& bi1 = boundary(SV1->info().id());
				Boundary& bi2 = boundary(SV2->info().id());

				Real d13 = bi1.p[bi1.coordinate] - (SV3->point())[bi1.coordinate];
				Real d23 = bi2.p[bi2.coordinate] - (SV3->point())[bi2.coordinate];
				cell->info().solidLine[j][facetF1] = abs(d23);
				cell->info().solidLine[j][facetF2] = abs(d13);
			}; break;
			default:
				// corner throat with 3 fictious vertices, no relevant solid line
				cell->info().solidLine[j][0] = cell->info().solidLine[j][1] = cell->info().solidLine[j][2] = 0;
		}

		lSolid = cell->info().solidLine[j][0] + cell->info().solidLine[j][1] + cell->info().solidLine[j][2];
		if (lSolid) cell->info().solidLine[j][3] = 1.0 / lSolid;
		else
			cell->info().solidLine[j][3] = 0;
	}
	template <class Tesselation> Real Network<Tesselation>::lineSolidFacet(Sphere ST1, Sphere ST2, Sphere ST3)
	{
		Real    line;
		Real    squaredRadius = ST1.weight();
		CVector v12 = ST2.point() - ST1.point();
		CVector v13 = ST3.point() - ST1.point();

		Real norme12 = v12.squared_length();
		Real norme13 = v13.squared_length();

		Real cosA = v12 * v13 / (sqrt(norme13 * norme12));
		line = acos(cosA) * sqrt(squaredRadius);
		return line;
	}

	template <class Tesselation> void Network<Tesselation>::setAlphaBoundary(Real alpha, bool fixed)
	{
		std::vector<Vector3r> vSegments;
		RTriangulation        temp(T[currentTes].Triangulation());
		AlphaShape            as(temp);
		alphaIdOffset = T[currentTes].vertexHandles.size(); // - 1;
		Real minAlpha = as.find_alpha_solid();

		if (alpha == 0) {
			alpha = minAlpha;
			as.set_alpha(alpha);
		} else {
			if (alpha < minAlpha) {
				cerr << "Using alpha<minAlpha will not work. Using minimum alpha" << endl;
				alpha = minAlpha;
			}
			as.set_alpha(alpha);
		}
		std::list<Facet> facets;
		as.get_alpha_shape_facets(std::back_inserter(facets), AlphaShape::REGULAR);
		for (auto fp = facets.begin(); fp != facets.end(); fp++) {
			Facet f = *fp;
			if (as.classify(f.first) != AlphaShape::INTERIOR) {
				f = as.mirror_facet(f);
				const CellHandle& outerCell = f.first->neighbor(f.second);
				if (as.is_infinite(outerCell) or fixed) {
					Sphere  sph;
					bool    b;
					CVector n;
					T[currentTes].circumCenter(f.first, f.second, alpha, b, sph, n);
					if (math::isnan(sph.point().x()) or math::isnan(sph.point().y()) or math::isnan(sph.point().z())) {
						cerr << "skipping isnan point" << endl;
						continue;
					}
					unsigned int id = T[currentTes].vertexHandles.size();
					VertexHandle Vh = T[currentTes].Triangulation().insert(Sphere(sph.point(), alpha)); //, pow(2*alphaRad-deltaAlpha,2)));
					if (Vh != NULL) {
						Vh->info().isAlpha = true;
						Vh->info().setId(id);
						Vh->info().isFictious = false;
						T[currentTes].vertexHandles.push_back(Vh);
						T[currentTes].maxId = math::max(T[currentTes].maxId, (int)id);
					} else {
						cerr << " : __Vh==NULL, convex__ :(" << endl;
					}
				} else {
					if (!outerCell->info().isAlpha) {
						outerCell->info().isAlpha = true;
						Point p = T[currentTes].setCircumCenter(outerCell);
						Real  weight
						        = (p - outerCell->vertex(0)->point().point()).squared_length() - outerCell->vertex(0)->point().weight();
						unsigned int id = T[currentTes].vertexHandles.size();
						VertexHandle Vh = T[currentTes].Triangulation().insert(Sphere(p, weight)); //pow(sqrt(weight)-deltaAlpha,2) ));
						if (Vh != NULL) {
							Vh->info().isAlpha = true;
							Vh->info().setId(id);
							Vh->info().isFictious = false;
						} else {
							cerr << " : __Vh==NULL, concave__ :(" << endl;
						}
						T[currentTes].vertexHandles.push_back(Vh);
						T[currentTes].maxId = math::max(T[currentTes].maxId, (int)id);
					}
				}
			}
		}
	}

} //namespace CGT

}; // namespace yade

#endif //FLOW_ENGINE
