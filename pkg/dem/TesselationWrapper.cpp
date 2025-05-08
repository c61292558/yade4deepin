/*************************************************************************
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_CGAL

#include "TesselationWrapper.hpp"
#include <lib/high-precision/Constants.hpp>
#include <preprocessing/dem/Shop.hpp>
#include <preprocessing/dem/SpherePack.hpp>

#ifdef YADE_OPENGL
#include <lib/opengl/GLUtils.hpp>
#include <lib/opengl/OpenGLWrapper.hpp>
#include <pkg/common/GLDrawFunctors.hpp>
#include <chrono>
#include <thread>
#endif

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((TesselationWrapper));
CREATE_LOGGER(TesselationWrapper);

// helper macro do assign Matrix3r values to subarrays
#define TENSOR_TO_MATRIX3R(mat, arr)                                                                                                                           \
	{                                                                                                                                                      \
		arr(0, 0) = mat(1, 1);                                                                                                                         \
		arr(0, 1) = mat(1, 2);                                                                                                                         \
		arr(0, 2) = mat(1, 3);                                                                                                                         \
		arr(1, 0) = mat(2, 1);                                                                                                                         \
		arr(1, 1) = mat(2, 2);                                                                                                                         \
		arr(1, 2) = mat(2, 3);                                                                                                                         \
		arr(2, 0) = mat(3, 1);                                                                                                                         \
		arr(2, 1) = mat(3, 2);                                                                                                                         \
		arr(2, 2) = mat(3, 3);                                                                                                                         \
	}

//spatial sort traits to use with a pair of CGAL::sphere pointers and integer.
//template<class _Triangulation>
struct RTraits_for_spatial_sort : public CGT::SimpleTriangulationTypes::RTriangulation::Geom_traits {
	typedef CGT::SimpleTriangulationTypes::RTriangulation::Geom_traits Gt;
	typedef std::pair<const CGT::Sphere*, Body::id_t>                  Point_3;

	struct Less_x_3 {
		bool operator()(const Point_3& p, const Point_3& q) const { return Gt::Less_x_3()(p.first->point(), q.first->point()); }
	};
	struct Less_y_3 {
		bool operator()(const Point_3& p, const Point_3& q) const { return Gt::Less_y_3()(p.first->point(), q.first->point()); }
	};
	struct Less_z_3 {
		bool operator()(const Point_3& p, const Point_3& q) const { return Gt::Less_z_3()(p.first->point(), q.first->point()); }
	};
	Less_x_3 less_x_3_object() const { return Less_x_3(); }
	Less_y_3 less_y_3_object() const { return Less_y_3(); }
	Less_z_3 less_z_3_object() const { return Less_z_3(); }
};


//function inserting points into a triangulation (where YADE::Sphere is converted to CGT::Sphere)
//and setting the info field to the bodies id.
//Possible improvements : use bodies pointers to avoid one copy, use aabb's lists to replace the shuffle/sort part
// template <class Triangulation>
void build_triangulation_with_ids(const shared_ptr<BodyContainer>& bodies, TesselationWrapper& TW, bool reset = true)
{
	if (reset) TW.clear();
	typedef SimpleTesselation::RTriangulation              RTriangulation;
	SimpleTesselation&                                     Tes = *(TW.Tes);
	RTriangulation&                                        T   = Tes.Triangulation();
	std::vector<CGT::Sphere>                               spheres;
	std::vector<std::pair<const CGT::Sphere*, Body::id_t>> pointsPtrs;
	spheres.reserve(bodies->size());
	pointsPtrs.reserve(bodies->size());
	Tes.vertexHandles.clear();
	Tes.vertexHandles.resize(bodies->size() + 6, NULL); //+6 extra slots in case boundaries will be added latter as additional vertices

	Body::id_t  Ng                = 0;
	Body::id_t& MaxId             = Tes.maxId;
	TW.mean_radius                = 0;
	int                nonSpheres = 0;
	shared_ptr<Sphere> sph(new Sphere);
	int                Sph_Index = sph->getClassIndexStatic();
	Scene*             scene     = Omega::instance().getScene().get();
	for (const auto& bi : *bodies) {
		if (bi and bi->shape->getClassIndex() == Sph_Index and bi->maskOk(TW.groupMask)) {
			const Sphere* s = YADE_CAST<Sphere*>(bi->shape.get());
			//FIXME: is the scene periodicity verification useful in the next line ? Tesselation seems to work in both periodic and non-periodic conditions with "scene->cell->wrapShearedPt(bi->state->pos)". I keep the verification to be consistent with all other uses of "wrapShearedPt" function.
			const Vector3r& pos = scene->isPeriodic ? scene->cell->wrapShearedPt(bi->state->pos) : bi->state->pos;
			const Real      rad = s->radius;
			CGT::Sphere     sp(CGT::Point(pos[0], pos[1], pos[2]), rad * rad);
			spheres.push_back(sp);
			pointsPtrs.push_back(std::make_pair(&(spheres[Ng] /*.point()*/), bi->getId()));
			TW.Pmin = CGT::Point(min(TW.Pmin.x(), pos.x() - rad), min(TW.Pmin.y(), pos.y() - rad), min(TW.Pmin.z(), pos.z() - rad));
			TW.Pmax = CGT::Point(max(TW.Pmax.x(), pos.x() + rad), max(TW.Pmax.y(), pos.y() + rad), max(TW.Pmax.z(), pos.z() + rad));
			Ng++;
			TW.mean_radius += rad;
			MaxId = max(MaxId, bi->getId());
		} else
			++nonSpheres;
	}
	TW.mean_radius /= Ng;
	TW.rad_divided = true;
	spheres.resize(Ng);
	pointsPtrs.resize(Ng);
	// random shuffling here is suggested in CGAL examples but it is probably not very helpful in yade since the positions are _usually_ random already.
	// We skip it to avoid unnecessary indeterminacy
	// std::random_shuffle(pointsPtrs.begin(), pointsPtrs.end());
	spatial_sort(pointsPtrs.begin(), pointsPtrs.end(), RTraits_for_spatial_sort() /*, CGT::RTriangulation::Weighted_point*/);

	RTriangulation::Cell_handle hint;

	TW.n_spheres = 0;
	for (std::vector<std::pair<const CGT::Sphere*, Body::id_t>>::const_iterator p = pointsPtrs.begin(); p != pointsPtrs.end(); ++p) {
		RTriangulation::Locate_type lt;
		RTriangulation::Cell_handle c;
		int                         li, lj;
		c                               = T.locate(*(p->first), lt, li, lj, hint);
		RTriangulation::Vertex_handle v = T.insert(*(p->first), lt, c, li, lj);
		if (v == RTriangulation::Vertex_handle()) hint = c;
		else {
			v->info().setId((unsigned int)p->second);
			//Vh->info().isFictious = false;//false is the default
			Tes.maxId                    = math::max(Tes.maxId, (int)p->second);
			Tes.vertexHandles[p->second] = v;
			hint                         = v->cell();
			++TW.n_spheres;
		}
	}
}

Real thickness = 0;

TesselationWrapper::~TesselationWrapper() { }

void TesselationWrapper::clear(void)
{
	Tes->Clear();
	Pmin        = CGT::Point(inf, inf, inf);
	Pmax        = CGT::Point(-inf, -inf, -inf);
	mean_radius = 0;
	n_spheres   = 0;
	rad_divided = false;
	bounded     = false;
	Tes->vertexHandles.clear();
}

void TesselationWrapper::insertSceneSpheres(bool reset) { build_triangulation_with_ids(scene->bodies, *this, reset); }

Real TesselationWrapper::Volume(unsigned int id) { return ((unsigned int)Tes->Max_id() >= id and Tes->vertex(id) != NULL) ? Tes->Volume(id) : 0; }

bool TesselationWrapper::insert(Real x, Real y, Real z, Real rad, unsigned int id)
{
	checkMinMax(x, y, z, rad);
	mean_radius += rad;
	++n_spheres;
	return (Tes->insert(x, y, z, rad, id) != NULL);
}

void TesselationWrapper::checkMinMax(Real x, Real y, Real z, Real rad)
{
	Pmin = CGT::Point(min(Pmin.x(), x - rad), min(Pmin.y(), y - rad), min(Pmin.z(), z - rad));
	Pmax = CGT::Point(max(Pmax.x(), x + rad), max(Pmax.y(), y + rad), max(Pmax.z(), z + rad));
}


bool TesselationWrapper::move(Real x, Real y, Real z, Real rad, unsigned int id)
{
	checkMinMax(x, y, z, rad);
	if (Tes->move(x, y, z, rad, id) != NULL) return true;
	else {
		std::cerr << "Tes->move(x,y,z,rad,id)==NULL" << std::endl;
		return false;
	}
}

void TesselationWrapper::computeTesselation(void)
{
	if (not(Tes->vertexHandles.size() > 0)) insertSceneSpheres();
	addBoundingPlanes();
	if (!rad_divided) {
		mean_radius /= n_spheres;
		rad_divided = true;
	}
	Tes->compute();
}

void TesselationWrapper::computeTesselation(Real pminx, Real pmaxx, Real pminy, Real pmaxy, Real pminz, Real pmaxz)
{
	if (not(Tes->vertexHandles.size() > 0)) insertSceneSpheres();
	addBoundingPlanes(pminx, pmaxx, pminy, pmaxy, pminz, pmaxz);
	computeTesselation();
}

void TesselationWrapper::computeVolumes(void)
{
	computeTesselation();
	Tes->computeVolumes();
}

int TesselationWrapper::addBoundingPlane(short axis, bool positive)
{
	Vector3r cornerMin   = Vector3r(Pmin.x(), Pmin.y(), Pmin.z());
	Vector3r cornerMax   = Vector3r(Pmax.x(), Pmax.y(), Pmax.z());
	Vector3r halfSize    = 0.5 * (cornerMax - cornerMin);
	Vector3r centerPoint = 0.5 * (cornerMin + cornerMax);
	// shift by half-size + a large radius
	Vector3r shift = Vector3r::Zero();
	shift[axis]    = positive ? 1 : -1;
	shift *= far * (cornerMax - cornerMin).norm();
	shift[axis] += positive ? halfSize[axis] : -halfSize[axis];
	centerPoint += shift;

	//find a free id
	int freeId = 0;
	while (Tes->vertexHandles[freeId] != NULL)
		++freeId;

	// we don't want to count this virtual sphere's radius in the average, so compute it before inserting
	if (!rad_divided) {
		mean_radius /= n_spheres;
		rad_divided = true;
	}
	// now insert
	Tes->vertexHandles[freeId] = Tes->insert(centerPoint[0], centerPoint[1], centerPoint[2], far * (cornerMax - cornerMin).norm(), freeId, false);
	return freeId;
}

void TesselationWrapper::addBoundingPlanes(Real pminx, Real pmaxx, Real pminy, Real pmaxy, Real pminz, Real pmaxz)
{
	if (!bounded) {
		if (!rad_divided) {
			mean_radius /= n_spheres;
			rad_divided = true;
		}
		// Insert the 6 additional vertices in the right place (usually they will be ids 0 to 5 when walls/facets/boxes are used, but not always)
		// append them at the end if the initial list is full
		int freeIds[6];
		int i = 0;
		for (int k = 0; k < 6; k++) {
			while (Tes->vertexHandles[i] != NULL)
				++i;
			freeIds[k] = i++;
		}
		// now insert
		Tes->vertexHandles[freeIds[0]] = Tes->insert(
		        pminx - far * (pmaxy - pminy), 0.5 * (pmaxy + pminy), 0.5 * (pmaxz + pminz), far * (pmaxy - pminy) + thickness, freeIds[0], true);
		Tes->vertexHandles[freeIds[1]] = Tes->insert(
		        pmaxx + far * (pmaxy - pminy), 0.5 * (pmaxy + pminy), 0.5 * (pmaxz + pminz), far * (pmaxy - pminy) + thickness, freeIds[1], true);
		Tes->vertexHandles[freeIds[2]] = Tes->insert(
		        0.5 * (pminx + pmaxx), pminy - far * (pmaxx - pminx), 0.5 * (pmaxz + pminz), far * (pmaxx - pminx) + thickness, freeIds[2], true);
		Tes->vertexHandles[freeIds[3]] = Tes->insert(
		        0.5 * (pminx + pmaxx), pmaxy + far * (pmaxx - pminx), 0.5 * (pmaxz + pminz), far * (pmaxx - pminx) + thickness, freeIds[3], true);
		Tes->vertexHandles[freeIds[4]] = Tes->insert(
		        0.5 * (pminx + pmaxx), 0.5 * (pmaxy + pminy), pminz - far * (pmaxy - pminy), far * (pmaxy - pminy) + thickness, freeIds[4], true);
		Tes->vertexHandles[freeIds[5]] = Tes->insert(
		        0.5 * (pminx + pmaxx), 0.5 * (pmaxy + pminy), pmaxz + far * (pmaxy - pminy), far * (pmaxy - pminy) + thickness, freeIds[5], true);

		bounded = true;
	}
}

void TesselationWrapper::addBoundingPlanes(void) { addBoundingPlanes(Pmin.x(), Pmax.x(), Pmin.y(), Pmax.y(), Pmin.z(), Pmax.z()); }

void TesselationWrapper::setState(bool state) { mma->setState(state ? 2 : 1, false, false, groupMask); }

void TesselationWrapper::loadState(string filename, bool stateNumber, bool bz2)
{
	CGT::TriaxialState& TS = stateNumber ? *(mma->analyser->TS1) : *(mma->analyser->TS0);
	TS.from_file(filename.c_str(), bz2);
}

void TesselationWrapper::saveState(string filename, bool stateNumber, bool bz2)
{
	CGT::TriaxialState& TS = stateNumber ? *(mma->analyser->TS1) : *(mma->analyser->TS0);
	TS.to_file(filename.c_str(), bz2);
}

void TesselationWrapper::defToVtkFromStates(string inputFile1, string inputFile2, string outputFile, bool bz2)
{
	mma->analyser->DefToFile(inputFile1.c_str(), inputFile2.c_str(), outputFile.c_str(), bz2);
}

void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool /*big*/, bool /*dynamic*/)
{
	body            = shared_ptr<Body>(new Body);
	body->groupMask = 2;
	shared_ptr<Sphere> iSphere(new Sphere);
	body->state->blockedDOFs = State::DOF_NONE;
	body->state->pos         = position;
	iSphere->radius          = radius;
	body->shape              = iSphere;
}

void TesselationWrapper::defToVtkFromPositions(string inputFile1, string inputFile2, string outputFile, bool /*bz2*/)
{
	SpherePack sp1, sp2;
	sp1.fromFile(inputFile1);
	sp2.fromFile(inputFile2);
	size_t imax = sp1.pack.size();
	if (imax != sp2.pack.size()) LOG_ERROR("The files have different numbers of spheres");
	shared_ptr<Body> body;
	scene->bodies->clear();
	for (size_t i = 0; i < imax; i++) {
		const SpherePack::Sph& sp(sp1.pack[i]);
		LOG_DEBUG("sphere (" << sp.c << " " << sp.r << ")");
		createSphere(body, sp.c, sp.r, false, true);
		scene->bodies->insert(body);
	}
	mma->setState(1);
	scene->bodies->clear();
	for (size_t i = 0; i < imax; i++) {
		const SpherePack::Sph& sp(sp2.pack[i]);
		LOG_DEBUG("sphere (" << sp.c << " " << sp.r << ")");
		createSphere(body, sp.c, sp.r, false, true);
		scene->bodies->insert(body);
	}
	mma->setState(2);
	mma->analyser->DefToFile(outputFile.c_str());
}

void TesselationWrapper::defToVtk(string outputFile) { mma->analyser->DefToFile(outputFile.c_str()); }

boost::python::dict TesselationWrapper::calcVolPoroDef(bool deformation)
{
	CGT::TriaxialState* ts;
	bounded = true; // TriaxialState already has bounding planes together with the actual spheres, no need to bound again.
	if (deformation) mma->analyser->computeParticlesDeformation();
	Tes = &mma->analyser->TS0->tesselation(); //no reason to use the final state if we don't want to compute deformations, keep using the initial
	ts  = mma->analyser->TS0;
	computeVolumes();
	RTriangulation& Tri = Tes->Triangulation();
	Pmin                = ts->box.base;
	Pmax                = ts->box.sommet;

	int              bodiesDim = Tes->Max_id() + 1; //=scene->bodies->size();
	vector<Real>     vol_(bodiesDim, 0);
	vector<Real>     poro_(bodiesDim, 0);
	vector<Matrix3r> def_(bodiesDim, Matrix3r::Zero());

	boost::python::list vol;
	boost::python::list poro;
	boost::python::list def;

	for (RTriangulation::Finite_vertices_iterator V_it = Tri.finite_vertices_begin(); V_it != Tri.finite_vertices_end(); V_it++) {
		const Body::id_t id = V_it->info().id();
		if (id < 0 or V_it->info().v() == 0 or V_it->info().isFictious) continue;
		Real sphereVol = 4.188790 * math::pow((V_it->point().weight()), 1.5); // 4/3*PI*R³ = 4.188...*R³
		vol_[id]       = V_it->info().v();
		poro_[id]      = (V_it->info().v() - sphereVol) / V_it->info().v();
		if (deformation) TENSOR_TO_MATRIX3R(mma->analyser->ParticleDeformation[id], def_[id]);
	}

	for (auto& v : vol_)
		vol.append(v);
	for (auto& v : poro_)
		poro.append(v);
	if (deformation)
		for (auto& v : def_)
			def.append(v);
	boost::python::dict ret;
	ret["vol"]  = vol;
	ret["poro"] = poro;
	if (deformation) ret["def"] = def;
	return ret;
}

boost::python::list TesselationWrapper::getAlphaFaces(Real alpha)
{
	vector<AlphaFace> faces;
	Tes->setAlphaFaces(faces, alpha);
	boost::python::list ret;
	for (auto f = faces.begin(); f != faces.end(); f++)
		ret.append(boost::python::make_tuple(Vector3i(f->ids[0], f->ids[1], f->ids[2]), makeVector3r(f->normal)));
	return ret;
}

boost::python::list TesselationWrapper::getAlphaCaps(Real alpha, Real shrinkedAlpha, bool fixedAlpha)
{
	vector<AlphaCap> caps;
	Tes->setExtendedAlphaCaps(caps, alpha, shrinkedAlpha, fixedAlpha);
	bounded = true;
	boost::python::list ret;
	for (auto f = caps.begin(); f != caps.end(); f++)
		ret.append(boost::python::make_tuple(f->id, makeVector3r(f->normal), makeVector3r(f->centroid)));
	return ret;
}

void TesselationWrapper::applyAlphaForces(Matrix3r stress, Real alpha, Real shrinkedAlpha, bool fixedAlpha, bool reset)
{
	if ((Tes->Triangulation().number_of_vertices() == 0) or reset) insertSceneSpheres(true);
	vector<AlphaCap> caps;
	Tes->setExtendedAlphaCaps(caps, alpha, shrinkedAlpha, fixedAlpha);
	bounded = true;
	for (const auto& b : *scene->bodies)
		scene->forces.setPermForce(b->id, Vector3r::Zero());
	for (auto f = caps.begin(); f != caps.end(); f++) {
		if (not scene->bodies->exists(f->id)) continue; // FIXME: probably not needed
		scene->forces.setPermForce(f->id, stress * makeVector3r(f->normal));
	}
}

void TesselationWrapper::applyAlphaVel(Matrix3r velGrad, Real alpha, Real shrinkedAlpha, bool fixedAlpha)
{
	build_triangulation_with_ids(scene->bodies, *this, true); //triangulation needed
	vector<AlphaCap> caps;
	Tes->setExtendedAlphaCaps(caps, alpha, shrinkedAlpha, fixedAlpha);
	bounded = true;
	for (const auto& b : *scene->bodies)
		b->state->blockedDOFs = State::DOF_NONE;
	const auto aabb     = Shop::aabbExtrema();
	Vector3r   bbCenter = 0.5 * (aabb.first + aabb.second);
	for (auto f = caps.begin(); f != caps.end(); f++) {
		Body* b               = Body::byId(f->id, scene).get();
		b->state->blockedDOFs = State::DOF_ALL;
		b->state->vel         = velGrad * (makeVector3r(f->centroid) - bbCenter);
		b->state->angVel      = Vector3r::Zero();
	}
}

Matrix3r TesselationWrapper::calcAlphaStress(Real alpha, Real shrinkedAlpha, bool fixedAlpha)
{
	build_triangulation_with_ids(scene->bodies, *this, true); //triangulation needed
	vector<AlphaCap> caps;
	Tes->setExtendedAlphaCaps(caps, alpha, shrinkedAlpha, fixedAlpha);
	bounded = true;
	Matrix3r cauchyLWS(Matrix3r::Zero());
	scene->forces.sync(); // needed to make resultants predictable
	alphaCapsVol = 0.;
	for (auto f = caps.begin(); f != caps.end(); f++) {
		Vector3r areaV     = makeVector3r(f->normal);
		Vector3r resultant = scene->forces.getPermForce(f->id) - scene->forces.getForce(f->id);
		Vector3r centroid  = makeVector3r(f->centroid);
		cauchyLWS += resultant * (centroid.transpose());
		alphaCapsVol += areaV.norm() / 3. * centroid.dot(areaV.normalized());
	}
	cauchyLWS /= alphaCapsVol;
	return cauchyLWS;
}

boost::python::list TesselationWrapper::getAlphaGraph(Real alpha, Real shrinkedAlpha, bool fixedAlpha)
{
	if (Tes->Triangulation().number_of_vertices() == 0) insertSceneSpheres(true);
	segments = Tes->getExtendedAlphaGraph(alpha, shrinkedAlpha, fixedAlpha);
	bounded  = true;
	boost::python::list ret;
	for (auto f = segments.begin(); f != segments.end(); f++)
		ret.append(*f);
	return ret;
}

boost::python::list TesselationWrapper::getAlphaVertices(Real alpha)
{
	vector<int>         vertices = Tes->getAlphaVertices(alpha);
	boost::python::list ret;
	for (auto f = vertices.begin(); f != vertices.end(); f++)
		ret.append(*f);
	return ret;
}

#ifdef YADE_OPENGL

YADE_PLUGIN((GlExtra_AlphaGraph))

GLUquadric* GlExtra_AlphaGraph::gluQuadric     = NULL;
int         GlExtra_AlphaGraph::glCylinderList = -1;
int         GlExtra_AlphaGraph::oneCylinder    = -1;

void GlExtra_AlphaGraph::render()
{
	if (not tesselationWrapper) tesselationWrapper = shared_ptr<TesselationWrapper>(new TesselationWrapper);
	if (tesselationWrapper->Tes->Triangulation().number_of_vertices() == 0) tesselationWrapper->insertSceneSpheres(true);
	if (tesselationWrapper->segments.size() == 0 or reset) {
#ifdef BREAK_OPENGL
		segments =
#endif
		        tesselationWrapper->segments = tesselationWrapper->Tes->getExtendedAlphaGraph(alpha, shrinkedAlpha, fixedAlpha);
		reset                                = true;
	}
#ifdef BREAK_OPENGL
	else
		segments = tesselationWrapper->segments;
#endif
	if ((reset or refreshDisplay) and not wire) {
		rots.clear();
		lengths.clear();
	}
	const vector<Vector3r>& segments_ = tesselationWrapper->segments;
	unsigned                maxI      = segments_.size() - 1;
	glLineWidth(lineWidth);
	glColor3v(color);
	if (lighting) glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
	glDisable(GL_CULL_FACE);
	if ((rots.size() == 0 or reset or refreshDisplay) and not wire) {
		Real     radiusTemp = 0;
		unsigned i          = 0;
		while (i < maxI) {
			const Vector3r& p1     = segments_[i];
			const Vector3r& p2     = segments_[i + 1];
			Vector3r        relPos = p2 - p1;
			Real            length = relPos.norm();
			lengths.push_back(length);
			rots.push_back(Eigen::Transform<Real, 3, Eigen::Affine>(Quaternionr().setFromTwoVectors(Vector3r(0, 0, 1), relPos / length)));
			radiusTemp += length;
			i += 2;
		}
		if (radius == 0) radius = radiusTemp / Real(i * 3);
		reset = false;
	}
	if (glCylinderList < 0) // will always be true unless the glGenLists() below is commented in (not sure it helps, should be tested)
	{
		unsigned i = 0, jj = 0;
		if (!gluQuadric) {
			gluQuadric = gluNewQuadric();
			if (!gluQuadric) throw runtime_error("Gl1_NormPhys::go unable to allocate new GLUquadric object (out of memory?).");
		}
		while (i < maxI) {
			const Vector3r& p1 = segments_[i];
			const Vector3r& p2 = segments_[i + 1];
			if (wire) {
				glBegin(GL_LINES)
					;
					glVertex3v(p1);
					glVertex3v(p2);
				glEnd();
			} else {
				glPushMatrix();
				glTranslatev(p1);
				glMultMatrix(rots[jj].data());
				gluCylinder(gluQuadric, radius, radius, lengths[jj], 4, 1);
				glPopMatrix();
			}
			i += 2;
			++jj;
		}
	}
}
#endif /*OPENGL*/

} // namespace yade

#endif /* YADE_CGAL */
