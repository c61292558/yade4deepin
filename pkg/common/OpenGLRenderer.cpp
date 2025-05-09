// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>

#ifdef YADE_OPENGL

#include "OpenGLRenderer.hpp"
#include <lib/high-precision/Constants.hpp>
#include <lib/opengl/GLUtils.hpp>
#include <lib/pyutil/gil.hpp>
#include <lib/serialization/EnumSupport.hpp>
#include <core/Aabb.hpp>
#include <core/Scene.hpp>
#include <core/Timing.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_ENUM(yade::OpenGLRenderer, BlinkHighlight, (NEVER)(NORMAL)(WEAK));
YADE_PLUGIN((OpenGLRenderer)(GlExtraDrawer));
CREATE_LOGGER(OpenGLRenderer);

void GlExtraDrawer::render()
{
	throw runtime_error("GlExtraDrawer::render called from class " + getClassName() + ". (did you forget to override it in the derived class?)");
}

bool      OpenGLRenderer::initDone = false;
const int OpenGLRenderer::numClipPlanes;
OpenGLRenderer::~OpenGLRenderer() { }

void OpenGLRenderer::init()
{
	for (const auto& item : Omega::instance().getDynlibsDescriptor()) {
		// if (Omega::instance().isInheritingFrom_recursive(item.first,"GlStateFunctor")) stateFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first, "GlBoundFunctor")) boundFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first, "GlShapeFunctor")) shapeFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first, "GlIGeomFunctor")) geomFunctorNames.push_back(item.first);
		if (Omega::instance().isInheritingFrom_recursive(item.first, "GlIPhysFunctor")) physFunctorNames.push_back(item.first);
	}
	initgl(); // creates functor objects in the proper sense

	clipPlaneNormals.resize(numClipPlanes);

	static bool glutInitDone = false;
	if (!glutInitDone) {
		glutInit(&Omega::instance().origArgc, Omega::instance().origArgv);
		/* transparent spheres (still not working): glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE | GLUT_ALPHA); glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE); */
		glutInitDone = true;
	}

	initDone = true;
}

void OpenGLRenderer::setBodiesRefSe3()
{
	LOG_DEBUG("(re)initializing reference positions and orientations.");
	const std::lock_guard<std::mutex> lock(scene->bodies->drawloopmutex);
	for (const auto& b : *scene->bodies)
		if (b && b->state) {
			b->state->refPos = b->state->pos;
			b->state->refOri = b->state->ori;
		}
	scene->cell->refHSize = scene->cell->hSize;
}

template <class FunctorType, class DispatcherT> void OpenGLRenderer::setupDispatcher(const vector<string>& names, DispatcherT& dispatcher)
{
	dispatcher.clearMatrix();
	for (const auto& s : names) {
		shared_ptr<FunctorType> f(boost::static_pointer_cast<FunctorType>(ClassFactory::instance().createShared(s)));
		f->initgl();
		dispatcher.add(f);
	}
}

void OpenGLRenderer::initgl()
{
	LOG_DEBUG("(re)initializing GL for gldraw methods.\n");
	setupDispatcher<GlBoundFunctor, GlBoundDispatcher>(boundFunctorNames, boundDispatcher);
	setupDispatcher<GlShapeFunctor, GlShapeDispatcher>(shapeFunctorNames, shapeDispatcher);
	setupDispatcher<GlIGeomFunctor, GlIGeomDispatcher>(geomFunctorNames, geomDispatcher);
	setupDispatcher<GlIPhysFunctor, GlIPhysDispatcher>(physFunctorNames, physDispatcher);
}

bool OpenGLRenderer::pointClipped(const Vector3r& p)
{
	if (numClipPlanes < 1) return false;
	for (int i = 0; i < numClipPlanes; i++)
		if (clipPlaneActive[i] && (p - clipPlaneSe3[i].position).dot(clipPlaneNormals[i]) < 0) return true;
	return false;
}

void OpenGLRenderer::setBodiesDispInfo()
{
	const std::lock_guard<std::mutex> lock(scene->bodies->drawloopmutex);
	if (scene->bodies->size() != bodyDisp.size()) {
		bodyDisp.resize(scene->bodies->size());
		for (unsigned k = 0; k < scene->bodies->size(); k++)
			bodyDisp[k].hidden = 0;
	}
	bool scaleRotations     = (rotScale != 1.0);
	bool scaleDisplacements = (dispScale != Vector3r::Ones());
	for (const auto& b : *scene->bodies) {
		if (!b || !b->state) continue;
		// ‘id’ shadows a member of ‘yade::OpenGLRenderer::id’
		size_t             id2     = b->getId();
		const Vector3r&    pos     = b->state->pos;
		const Vector3r&    refPos  = b->state->refPos;
		const Quaternionr& ori     = b->state->ori;
		const Quaternionr& refOri  = b->state->refOri;
		Vector3r           cellPos = (!scene->isPeriodic ? pos : scene->cell->wrapShearedPt(pos)); // inside the cell if periodic, same as pos otherwise
		bodyDisp[id2].isDisplayed  = !pointClipped(cellPos);
		// if no scaling and no periodic, return quickly
		if (!(scaleDisplacements || scaleRotations || scene->isPeriodic)) {
			bodyDisp[id2].pos = pos;
			bodyDisp[id2].ori = ori;
			continue;
		}
		// apply scaling
		bodyDisp[id2].pos = cellPos; // point of reference (inside the cell for periodic)
		if (scaleDisplacements) bodyDisp[id2].pos += dispScale.cwiseProduct(Vector3r(pos - refPos)); // add scaled translation to the point of reference
		if (!scaleRotations) bodyDisp[id2].ori = ori;
		else {
			Quaternionr relRot = refOri.conjugate() * ori;
			AngleAxisr  aa(relRot);
			aa.angle() *= rotScale;
			bodyDisp[id2].ori = refOri * Quaternionr(aa);
		}
	}
}

// draw periodic cell, if active
void OpenGLRenderer::drawPeriodicCell()
{
	if (!scene->isPeriodic) return;
	glColor3v(cellColor);
	glPushMatrix();
	// Vector3r size=scene->cell->getSize();
	const Matrix3r& hSize = scene->cell->hSize;
	if (dispScale != Vector3r::Ones()) {
		const Matrix3r& refHSize(scene->cell->refHSize);
		Matrix3r        scaledHSize;
		for (int i = 0; i < 3; i++)
			scaledHSize.col(i) = refHSize.col(i) + dispScale.cwiseProduct(Vector3r(hSize.col(i) - refHSize.col(i)));
		GLUtils::Parallelepiped(scaledHSize.col(0), scaledHSize.col(1), scaledHSize.col(2));
	} else {
		GLUtils::Parallelepiped(hSize.col(0), hSize.col(1), hSize.col(2));
	}
	glPopMatrix();
}

void OpenGLRenderer::resetSpecularEmission()
{
	glMateriali(GL_FRONT, GL_SHININESS, 80);
	const GLfloat glutMatSpecular[4] = { 0.3f, 0.3f, 0.3f, 0.5f };
	const GLfloat glutMatEmit[4]     = { 0.2f, 0.2f, 0.2f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, glutMatSpecular);
	glMaterialfv(GL_FRONT, GL_EMISSION, glutMatEmit);
}

void OpenGLRenderer::render(const shared_ptr<Scene>& _scene, Body::id_t selection)
{
	gilLock lockgil;
	if (!initDone) init();
	assert(initDone);
	selId = selection;

	scene = _scene;

	// assign scene inside functors
	boundDispatcher.updateScenePtr();
	geomDispatcher.updateScenePtr();
	physDispatcher.updateScenePtr();
	shapeDispatcher.updateScenePtr();
	// stateDispatcher.updateScenePtr();

	// just to make sure, since it is not initialized by default
	if (!scene->bound) scene->bound = shared_ptr<Aabb>(new Aabb);

	// recompute emissive light colors for highlighted bodies
	Real now              = TimingInfo::getNow(/*even if timing is disabled*/ true) * 1e-9;
	highlightEmission0[0] = highlightEmission0[1] = highlightEmission0[2] = ((blinkHighlight == BlinkHighlight::WEAK) ? 0.2 : 0.8) * normSquare(now, 1);
	highlightEmission1[0] = highlightEmission1[1] = highlightEmission1[2] = ((blinkHighlight == BlinkHighlight::WEAK) ? 0.2 : 0.6) * normSaw(now, 2);

	// clipping planes
	assert(clipPlaneNormals.size() == (size_t)numClipPlanes);
	for (size_t i = 0; i < numClipPlanes; ++i) {
		// someone could have modified those from python and truncate the vectors; fill those here in that case
		if (i == clipPlaneSe3.size()) clipPlaneSe3.emplace_back(Se3r(Vector3r::Zero(), Quaternionr::Identity()));
		if (i == clipPlaneActive.size()) clipPlaneActive.emplace_back(0);
		if (i == clipPlaneNormals.size()) clipPlaneNormals.emplace_back(Vector3r::UnitX());

		// end filling stuff modified from python
		if (clipPlaneActive[i]) clipPlaneNormals[i] = clipPlaneSe3[i].orientation * Vector3r(0, 0, 1);
	}
	// set displayed Se3 of body (scaling) and isDisplayed (clipping)
	setBodiesDispInfo();

	glClearColor(bgColor[0], bgColor[1], bgColor[2], 1.0);

	// set light sources
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1); // important: do lighting calculations on both sides of polygons

	const GLfloat pos[4]           = { (float)lightPos[0], (float)lightPos[1], (float)lightPos[2], 1.0 };
	const GLfloat ambientColor[4]  = { 0.2f, 0.2f, 0.2f, 1.0f };
	const GLfloat specularColor[4] = { 1, 1, 1, 1.f };
	const GLfloat diffuseLight[4]  = { (float)lightColor[0], (float)lightColor[1], (float)lightColor[2], 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	if (light1) glEnable(GL_LIGHT0);
	else
		glDisable(GL_LIGHT0);

	const GLfloat pos2[4]           = { (float)light2Pos[0], (float)light2Pos[1], (float)light2Pos[2], 1.0 };
	const GLfloat ambientColor2[4]  = { 0.0, 0.0, 0.0, 1.0 };
	const GLfloat specularColor2[4] = { 1.f, 1.f, 0.6f, 1.f };
	const GLfloat diffuseLight2[4]  = { (float)light2Color[0], (float)light2Color[1], (float)light2Color[2], 1.0f };
	glLightfv(GL_LIGHT1, GL_POSITION, pos2);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specularColor2);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambientColor2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight2);
	if (light2) glEnable(GL_LIGHT1);
	else
		glDisable(GL_LIGHT1);

	glEnable(GL_LIGHTING);

	glEnable(GL_CULL_FACE);
	// http://www.sjbaker.org/steve/omniv/opengl_lighting.html
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	//Shared material settings
	resetSpecularEmission();

	drawPeriodicCell();

	if (dof || id) renderDOF_ID();
	if (bound) renderBound();

	if (shape) renderShape();
	if (intrAllWire) renderAllInteractionsWire();
	if (intrGeom) renderIGeom();
	if (intrPhys) renderIPhys();

	for (const auto& d : extraDrawers) {
		if (d->dead) continue;
		glPushMatrix();
		d->scene = scene.get();
		d->render();
		glPopMatrix();
	}
}

void OpenGLRenderer::renderAllInteractionsWire()
{
	const std::lock_guard<std::mutex> lockB(scene->bodies->drawloopmutex);
	const std::lock_guard<std::mutex> lock(scene->interactions->drawloopmutex);
	for (const auto& i : *scene->interactions) {
		// geometry must exist              , sometimes a body can get deleted
		if ((not i->functorCache.geomExists)) { continue; }
		const boost::shared_ptr<const Body> b1 = Body::byId(i->getId1(), scene);
		const boost::shared_ptr<const Body> b2 = Body::byId(i->getId2(), scene);
		// If the Body gets deleted after the two lines above, then we hold the last instance of it. And it will be deleted when this part of the loop ends
		// using `const` makes this threadsafe, because shared_ptr holds a threadsafe atomic count of instances, while const means that we are not writing there.
		// Only reading the soon-to-be-deleted (in next 50 miliseconds) body position.
		if ((not b1) or (not b2)) { continue; }
		glColor3v(i->isReal() ? Vector3r(0, 1, 0) : Vector3r(.5, 0, 1));
		Vector3r        p1   = b1->state->pos;
		const Vector3r& size = scene->cell->getSize();
		Vector3r        shift2(i->cellDist[0] * size[0], i->cellDist[1] * size[1], i->cellDist[2] * size[2]);
		// in sheared cell, apply shear on the mutual position as well
		shift2       = scene->cell->shearPt(shift2);
		Vector3r rel = b2->state->pos + shift2 - p1;
		if (scene->isPeriodic) p1 = scene->cell->wrapShearedPt(p1);
		glBegin(GL_LINES)
			;
			glVertex3v(p1);
			glVertex3v(Vector3r(p1 + rel));
		glEnd();
	}
}

void OpenGLRenderer::renderDOF_ID()
{
	const GLfloat                     ambientColorSelected[4]   = { 10.0, 0.0, 0.0, 1.0 };
	const GLfloat                     ambientColorUnselected[4] = { 0.5, 0.5, 0.5, 1.0 };
	const std::lock_guard<std::mutex> lock(scene->bodies->drawloopmutex);
	for (const auto& b : *scene->bodies) {
		if (!b) continue;
		if (b->shape && ((b->getGroupMask() & mask) || b->getGroupMask() == 0)) {
			if (!id && b->state->blockedDOFs == 0) continue;
			if (selId == b->getId()) { glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColorSelected); }
			{ // write text
				glColor3(1.0 - bgColor[0], 1.0 - bgColor[1], 1.0 - bgColor[2]);
				unsigned    d    = b->state->blockedDOFs;
				std::string sDof = std::string() + (((d & State::DOF_X) != 0) ? "x" : "") + (((d & State::DOF_Y) != 0) ? "y" : " ")
				        + (((d & State::DOF_Z) != 0) ? "z" : "") + (((d & State::DOF_RX) != 0) ? "X" : "")
				        + (((d & State::DOF_RY) != 0) ? "Y" : "") + (((d & State::DOF_RZ) != 0) ? "Z" : "");
				std::string sId = boost::lexical_cast<std::string>(b->getId());
				std::string str;
				if (dof && id) sId += " ";
				if (id) str += sId;
				if (dof) str += sDof;
				const Vector3r& h(selId == b->getId() ? highlightEmission0 : Vector3r(1, 1, 1));
				glColor3v(h);
				GLUtils::GLDrawText(str, bodyDisp[b->id].pos, h);
			}
			if (selId == b->getId()) { glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColorUnselected); }
		}
	}
}

void OpenGLRenderer::renderIGeom()
{
	geomDispatcher.scene = scene.get();
	geomDispatcher.updateScenePtr();
	{
		const std::lock_guard<std::mutex> lockB(scene->bodies->drawloopmutex);
		const std::lock_guard<std::mutex> lock(scene->interactions->drawloopmutex);
		for (const auto& I : *scene->interactions) {
			if (!I->geom) continue;        // avoid refcount manipulations if the interaction is not real anyway
			shared_ptr<IGeom> ig(I->geom); // keep reference so that ig does not disappear suddenly while being rendered
			if (!ig) continue;
			const auto b1 = Body::byId(I->getId1(), scene), b2 = Body::byId(I->getId2(), scene);
			if (!(bodyDisp[I->getId1()].isDisplayed || bodyDisp[I->getId2()].isDisplayed)) continue;
			glPushMatrix();
			geomDispatcher(ig, I, b1, b2, intrWire);
			glPopMatrix();
		}
	}
}

void OpenGLRenderer::renderIPhys()
{
	physDispatcher.scene = scene.get();
	physDispatcher.updateScenePtr();
	{
		const std::lock_guard<std::mutex> lockB(scene->bodies->drawloopmutex);
		const std::lock_guard<std::mutex> lock(scene->interactions->drawloopmutex);
		for (const auto& I : *scene->interactions) {
			shared_ptr<IPhys> ip(I->phys);
			if (!ip) continue;
			const auto b1 = Body::byId(I->getId1(), scene), b2 = Body::byId(I->getId2(), scene);
			Body::id_t id1 = I->getId1(), id2 = I->getId2();
			if (!(bodyDisp[id1].isDisplayed || bodyDisp[id2].isDisplayed)) continue;
			glPushMatrix();
			physDispatcher(ip, I, b1, b2, intrWire);
			glPopMatrix();
		}
	}
}

void OpenGLRenderer::renderBound()
{
	boundDispatcher.scene = scene.get();
	boundDispatcher.updateScenePtr();

	const std::lock_guard<std::mutex> lock(scene->bodies->drawloopmutex);
	for (const auto& b : *scene->bodies) {
		if (!b || !b->bound) continue;
		if (!bodyDisp[b->getId()].isDisplayed or bodyDisp[b->getId()].hidden) continue;
		if (b->bound && ((b->getGroupMask() & mask) || b->getGroupMask() == 0)) {
			glPushMatrix();
			boundDispatcher(b->bound, scene.get());
			glPopMatrix();
		}
	}
	// since we remove the functor as Scene doesn't inherit from Body anymore, hardcore the rendering routine here
	// for periodic scene, renderPeriodicCell is called separately
	if (!scene->isPeriodic) {
		if (!scene->bound) scene->updateBound();
		glColor3v(Vector3r(0, 1, 0));
		Vector3r size   = scene->bound->max - scene->bound->min;
		Vector3r center = .5 * (scene->bound->min + scene->bound->max);
		glPushMatrix();
		glTranslatev(center);
		glScalev(size);
		glutWireCube(1);
		glPopMatrix();
	}
}

// this function is called for both rendering as well as
// in the selection mode

// nice reading on OpenGL selection
// http://glprogramming.com/red/chapter13.html

void OpenGLRenderer::renderShape()
{
	shapeDispatcher.scene = scene.get();
	shapeDispatcher.updateScenePtr();

	// instead of const shared_ptr&, get proper shared_ptr;
	// Less efficient in terms of performance, since memory has to be written (not measured, though),
	// but it is still better than crashes if the body gets deleted meanwile.
	const std::lock_guard<std::mutex> lock(scene->bodies->drawloopmutex);
	for (const auto& b : *scene->bodies) {
		if (!b || !b->shape) continue;
		if (!(bodyDisp[b->getId()].isDisplayed and !bodyDisp[b->getId()].hidden)) continue;
		Vector3r    pos = bodyDisp[b->getId()].pos;
		Quaternionr ori = bodyDisp[b->getId()].ori;
		if (!b->shape || !((b->getGroupMask() & mask) || b->getGroupMask() == 0)) continue;

		// ignored in non-selection mode, use it always
		glPushName(b->id);
		bool highlight
		        = (b->id == selId || (b->clumpId >= 0 && b->clumpId == selId) || b->shape->highlight) and (blinkHighlight != BlinkHighlight::NEVER);

		glPushMatrix();
		AngleAxisr aa(ori);
		glTranslate(pos[0], pos[1], pos[2]);
		glRotate(aa.angle() * Mathr::RAD_TO_DEG, aa.axis()[0], aa.axis()[1], aa.axis()[2]);
		if (highlight) {
			// set hightlight
			// different color for body highlighted by selection and by the shape attribute
			const Vector3r& h((selId == b->id || (b->clumpId >= 0 && selId == b->clumpId)) ? highlightEmission0 : highlightEmission1);
			glMaterialv(GL_FRONT_AND_BACK, GL_EMISSION, h);
			glMaterialv(GL_FRONT_AND_BACK, GL_SPECULAR, h);
			shapeDispatcher(b->shape, b->state, wire || b->shape->wire, viewInfo);
			// reset highlight
			resetSpecularEmission();
		} else {
			// no highlight; in case previous functor fiddled with glMaterial
			resetSpecularEmission();
			shapeDispatcher(b->shape, b->state, wire || b->shape->wire, viewInfo);
		}
		glPopMatrix();
		if (highlight) {
			if (!b->bound || wire || b->shape->wire) GLUtils::GLDrawInt(b->getId(), pos);
			else {
				// move the label towards the camera by the bounding box so that it is not hidden inside the body
				const Vector3r& mn = b->bound->min;
				const Vector3r& mx = b->bound->max;
				const Vector3r& p  = pos;
				Vector3r        ext(
                                        viewDirection[0] > 0 ? p[0] - mn[0] : p[0] - mx[0],
				        viewDirection[1] > 0 ? p[1] - mn[1] : p[1] - mx[1],
				        viewDirection[2] > 0 ? p[2] - mn[2] : p[2] - mx[2]); // signed extents towards the camera
				Vector3r dr = -1.01 * (viewDirection.dot(ext) * viewDirection);
				GLUtils::GLDrawInt(b->getId(), pos + dr, Vector3r::Ones());
			}
		}
		// if the body goes over the cell margin, draw it in positions where the bbox overlaps with the cell in wire
		// precondition: pos is inside the cell.
		if (b->bound && scene->isPeriodic && ghosts) {
			const Vector3r& cellSize(scene->cell->getSize());
			pos = scene->cell->unshearPt(pos); // remove the shear component
			// traverse all periodic cells around the body, to see if any of them touches
			Vector3r halfSize = b->bound->max - b->bound->min;
			halfSize *= .5;
			Vector3r pmin, pmax;
			Vector3i i;
			for (i[0] = -1; i[0] <= 1; i[0]++)
				for (i[1] = -1; i[1] <= 1; i[1]++)
					for (i[2] = -1; i[2] <= 1; i[2]++) {
						if (i[0] == 0 && i[1] == 0 && i[2] == 0) continue; // middle; already rendered above
						Vector3r pos2 = pos
						        + Vector3r(cellSize[0] * i[0], cellSize[1] * i[1], cellSize[2] * i[2]); // shift, but without shear!
						pmin = pos2 - halfSize;
						pmax = pos2 + halfSize;
						if (pmin[0] <= cellSize[0] && pmax[0] >= 0 && pmin[1] <= cellSize[1] && pmax[1] >= 0 && pmin[2] <= cellSize[2]
						    && pmax[2] >= 0) {
							Vector3r pt = scene->cell->shearPt(pos2);
							if (pointClipped(pt)) continue;
							glLoadName(b->id);
							glPushMatrix();
							glTranslatev(pt);
							glRotate(aa.angle() * Mathr::RAD_TO_DEG, aa.axis()[0], aa.axis()[1], aa.axis()[2]);
							shapeDispatcher(b->shape, b->state, /*wire*/ true, viewInfo);
							glPopMatrix();
						}
					}
		}
		glPopName();
	}
}

} // namespace yade

#endif /* YADE_OPENGL */
