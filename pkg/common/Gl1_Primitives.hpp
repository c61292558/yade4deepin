/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  Copyright (C) 2008 by Sergei Dorofeenko                               *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  © 2008 Václav Šmilauer <eudoxos@arcig.cz>                             *
*                                                                        *
*  © 2008 Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <core/Aabb.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/GLDrawFunctors.hpp>
#include <pkg/common/Sphere.hpp>

namespace yade { // Cannot have #include directive inside.

class Gl1_Aabb : public GlBoundFunctor {
public:
	void go(const shared_ptr<Bound>&, Scene*) override;
	RENDERS(Aabb);
	// clang-format off
	YADE_CLASS_BASE_DOC(Gl1_Aabb,GlBoundFunctor,"Render Axis-aligned bounding box (:yref:`Aabb`).");
	// clang-format on
};
REGISTER_SERIALIZABLE(Gl1_Aabb);


class Gl1_Box : public GlShapeFunctor {
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	RENDERS(Box);
	// clang-format off
	YADE_CLASS_BASE_DOC(Gl1_Box,GlShapeFunctor,"Renders :yref:`Box` object");
	// clang-format on
};

REGISTER_SERIALIZABLE(Gl1_Box);

class Gl1_Facet : public GlShapeFunctor {
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	RENDERS(Facet);
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_Facet,GlShapeFunctor,"Renders :yref:`Facet` object",
		((bool,normals,false,,"In wire mode, render normals of facets and edges; facet's :yref:`colors<Shape::color>` are disregarded in that case."))
	);
	// clang-format on
};

REGISTER_SERIALIZABLE(Gl1_Facet);


class Gl1_Sphere : public GlShapeFunctor {
private:
	// for stripes
	static vector<Vector3r> vertices, faces;
	static int              glStripedSphereList;
	static int              glGlutSphereList;
	void                    subdivideTriangle(Vector3r& v1, Vector3r& v2, Vector3r& v3, int depth);
	//Generate GlList for GLUT sphere
	void initGlutGlList();
	//Generate GlList for sliced spheres
	void initStripedGlList();
	//for regenerating glutSphere or glutTorus list if needed
	static Real prevQuality;
	//for regenerating glutSphere or glutTorus list if needed
	static string prevDisplayMode;
	//for regenerating glutTorus list if needed
	static char prevCircleAllowedRotationAxis;

public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_Sphere,GlShapeFunctor,"Renders :yref:`Sphere` object",
		((Real,quality,1.0,,"Change discretization level of spheres. quality>1  for better image quality, at the price of more cpu/gpu usage, 0<quality<1 for faster rendering. If mono-color spheres are displayed (:yref:`Gl1_Sphere::stripes` = False), quality mutiplies :yref:`Gl1_Sphere::glutSlices` and :yref:`Gl1_Sphere::glutStacks`. If striped spheres are displayed (:yref:`Gl1_Sphere::stripes` = True), only integer increments are meaningfull : quality=1 and quality=1.9 will give the same result, quality=2 will give finer result."))
		((bool,wire,false,,"Only show wireframe (controlled by ``glutSlices`` and ``glutStacks``."))
		((bool,stripes,true,,"In non-wire rendering, show stripes clearly showing particle rotation."))
		((bool,localSpecView,true,,"Compute specular light in local eye coordinate system."))
		((int,glutSlices,12,(Attr::noSave | Attr::readonly),"Base number of sphere slices, multiplied by :yref:`Gl1_Sphere::quality` before use); not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
		((int,glutStacks,6,(Attr::noSave | Attr::readonly),"Base number of sphere stacks, multiplied by :yref:`Gl1_Sphere::quality` before use; not used with ``stripes`` (see `glut{Solid,Wire}Sphere reference <http://www.opengl.org/documentation/specs/glut/spec3/node81.html>`_)"))
		((bool,circleView,false,,"For 2D simulations : display tori instead of spheres, so they will appear like circles if the viewer is looking in the right direction. In this case, remember to disable perspective by pressing \"t\"-key in the viewer."))
		((Real,circleRelThickness,0.2,,"If :yref:`Gl1_Sphere::circleView` is enabled, this is the torus diameter relative to the sphere radius (i.e. the circle relative thickness)."))
		((char,circleAllowedRotationAxis,'z',,"If :yref:`Gl1_Sphere::circleView` is enabled, this is the only axis ('x', 'y' or 'z') along which rotation is allowed for the 2D simulation. It allows right orientation of the tori to appear like circles in the viewer. For example, if circleAllowedRotationAxis='x' is set, blockedDOFs=\"YZ\" should also be set for all your particles."))
	);
	// clang-format on
	RENDERS(Sphere);
};

REGISTER_SERIALIZABLE(Gl1_Sphere);

} // namespace yade
