/*************************************************************************
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef YADE_CGAL

#pragma once
#include <lib/triangulation/Tesselation.h>
#include <core/GlobalEngine.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/MicroMacroAnalyser.hpp>
#ifdef YADE_OPENGL
#include <lib/opengl/GLUtils.hpp>
#include <lib/opengl/OpenGLWrapper.hpp>
#include <pkg/common/OpenGLRenderer.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

/*! \class TesselationWrapper
 * \brief Handle the triangulation of spheres in a scene, build tesselation on request, and give access to computed quantities : currently volume and porosity of each Voronoï sphere.
 * More accessors in course of implementation. Feel free to suggest new ones.
 *
 * Example usage script :
 *
 *tt=TriaxialTest()
 *tt.generate("test.yade")
 *O.load("test.yade")
 *O.run(100,True) 
 *TW=TesselationWrapper()
 *TW.triangulate() #compute regular Delaunay triangulation, don't construct tesselation
 *TW.computeVolumes() #will silently tesselate the packing
 *TW.volume(10) #get volume associated to sphere of id 10
 *
 */


class TesselationWrapper : public GlobalEngine {
public:
	typedef CGT::_Tesselation<CGT::SimpleTriangulationTypes> Tesselation;
	typedef Tesselation::RTriangulation                      RTriangulation;
	typedef Tesselation::VertexInfo                          VertexInfo;
	typedef Tesselation::CellInfo                            CellInfo;
	typedef RTriangulation::Finite_edges_iterator            FiniteEdgesIterator;
	typedef Tesselation::AlphaFace                           AlphaFace;
	typedef Tesselation::AlphaCap                            AlphaCap;

	Tesselation      tesObj;
	Tesselation*     Tes;
	Real             mean_radius, inf;
	bool             rad_divided;
	bool             bounded;
	CGT::Point       Pmin;
	CGT::Point       Pmax;
	vector<Vector3r> segments;

	~TesselationWrapper();

	/// Insert a sphere, "id" will be used by some getters to retrieve spheres
	bool insert(Real x, Real y, Real z, Real rad, unsigned int id);
	/// A faster version of insert, inserting all spheres in scene (first erasing current triangulation  if reset=true)
	void insertSceneSpheres(bool reset = true);
	/// Move one sphere to the new position (x,y,z) and maintain triangulation (invalidates the tesselation)
	bool move(Real x, Real y, Real z, Real rad, unsigned int id);

	void checkMinMax(Real x, Real y, Real z, Real rad); //for experimentation purpose
	                                                    /// Reset the triangulation
	void clear(void);

	/// Add 6 axis-aligned bounding planes (modelised as spheres with (almost) infinite radius)
	void addBoundingPlanes(void);
	/// Force boudaries at positions not equal to precomputed ones
	void addBoundingPlanes(Real pminx, Real pmaxx, Real pminy, Real pmaxy, Real pminz, Real pmaxz);
	/// Insert one single axis-aligned bounding plane
	int addBoundingPlane(short axis, bool positive);

	///compute voronoi centers then stop (don't compute anything else)
	void computeTesselation(void);
	void computeTesselation(Real pminx, Real pmaxx, Real pminy, Real pmaxy, Real pminz, Real pmaxz);

	void                testAlphaShape(Real alpha) { Tes->testAlphaShape(alpha); }
	boost::python::list getAlphaFaces(Real alpha);
	boost::python::list getAlphaCaps(Real alpha, Real shrinkedAlpha, bool fixedAlpha);
	boost::python::list getAlphaVertices(Real alpha);
	boost::python::list getAlphaGraph(Real alpha, Real shrinkedAlpha, bool fixedAlpha);
	void                applyAlphaForces(Matrix3r stress, Real alpha, Real shrinkedAlpha, bool fixedAlpha, bool reset = true);
	void                applyAlphaVel(Matrix3r velGrad, Real alpha, Real shrinkedAlpha, bool fixedAlpha);
	Matrix3r            calcAlphaStress(Real alpha, Real shrinkedAlpha, bool fixedAlpha);

	///compute Voronoi vertices + volumes of all cells
	///use computeTesselation to force update, e.g. after spheres positions have been updated
	void computeVolumes(void);
	void computeDeformations(void) { mma->analyser->computeParticlesDeformation(); }
	///Get volume of the sphere inserted with indentifier "id""
	Real Volume(unsigned int id);
	Real deformation(unsigned int id, unsigned int i, unsigned int j)
	{
		if (!mma->analyser->ParticleDeformation.size()) {
			LOG_ERROR("compute deformations first");
			return 0;
		}
		if (mma->analyser->ParticleDeformation.size() < id) {
			LOG_ERROR("id out of bounds");
			return 0;
		}
		if (i < 1 || i > 3 || j < 1 || j > 3) {
			LOG_ERROR("tensor index must be between 1 and 3");
			return 0;
		}
		return mma->analyser->ParticleDeformation[id](i, j);
	}

	Matrix3r deformationTensor(unsigned int id)
	{
		if (!mma->analyser->ParticleDeformation.size()) {
			LOG_ERROR("compute deformations first");
			return Matrix3r::Zero();
		}
		if (mma->analyser->ParticleDeformation.size() < id) {
			LOG_ERROR("id out of bounds");
			return Matrix3r::Zero();
		}
		auto     m = mma->analyser->ParticleDeformation[id];
		Matrix3r m3;
		m3 << m(1, 1), m(1, 2), m(1, 3), m(2, 1), m(2, 2), m(2, 3), m(3, 1), m(3, 2), m(3, 3);
		return m3;
	}

	/// make the current state the initial (0) or final (1) configuration for the definition of displacement increments, use only state=0 if you just want to get only volmumes and porosity
	void setState(bool state = 0);
	void loadState(string fileName, bool stateNumber = 0, bool bz2 = false);
	void saveState(string fileName, bool stateNumber = 0, bool bz2 = false);
	/// read two state files and write per-particle deformation to a vtk file. The second variant uses existing states.
	void defToVtkFromStates(string inputFile1, string inputFile2, string outputFile = "def.vtk", bool bz2 = false);
	void defToVtkFromPositions(string inputFile1, string inputFile2, string outputFile = "def.vtk", bool bz2 = false);
	void defToVtk(string outputFile = "def.vtk");

	/// return python array containing voronoi volumes, per-particle porosity, and optionaly per-particle deformation, if states 0 and 1 have been assigned
	boost::python::dict calcVolPoroDef(bool deformation);


public:
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(TesselationWrapper,GlobalEngine,"Handle the triangulation of spheres in a scene, build tesselation on request, and give access to computed quantities (see also the :ref:`dedicated section in user manual <MicroStressAndMicroStrain>`). The calculation of microstrain is explained in [Catalano2014a]_ \n\nSee example usage in script example/tesselationWrapper/tesselationWrapper.py.\n\nBelow is an output of the :yref:`defToVtk<TesselationWrapper::defToVtk>` function visualized with paraview (in this case Yade's TesselationWrapper was used to process experimental data obtained on sand by Edward Ando at Grenoble University, 3SR lab.)\n\n.. figure:: fig/localstrain.*\n\t:width: 9cm\n\nThe definition of outer contours of arbitrary shapes and the application of stress on them, based on CGAL's 'alpha shapes' is also possible. See :ysrc:`scripts/examples/alphaShapes/GlDrawAlpha.py` (giving the figure below) and other examples therein. Read more in [Pekmezi2020]_ and further papers by the same authors. \n\n.. figure:: fig/alphaShape.*\n\t:width: 9cm",
	((unsigned int,n_spheres,0,,"|ycomp|"))
	((Real,far,10000.,,"Defines the radius of the large virtual spheres used to define nearly flat boundaries around the assembly. The radius will be the (scene's) bounding box size multiplied by 'far'. Higher values will minimize the error theoretically (since the infinite sphere really defines a plane), but it may increase numerical errors at some point. The default should give a resonable compromize."))
	((Real,alphaCapsVol,0.,,"The volume of the packing as defined by the boundary alpha cap polygons"))
	((Matrix3r,grad_u,Matrix3r::Zero(),,"The Displacement Gradient Tensor"))
	((mask_t,groupMask,0,,"Bitmask for filtering spheres, ignored if 0."))
	((shared_ptr<MicroMacroAnalyser>, mma, new MicroMacroAnalyser,, "underlying object processing the data - see specific settings in :yref:`MicroMacroAnalyser` class documentation"))
	,/*deprec*/
	,/*init*/
	,/*ctor*/
	Tes = &tesObj;
	inf=1e10;
	mma->analyser->SetConsecutive(false);
	,/*py*/
	.def("triangulate",&TesselationWrapper::insertSceneSpheres,(boost::python::arg("reset")=true),"triangulate spheres of the packing")
	.def("addBoundingPlane",&TesselationWrapper::addBoundingPlane,((boost::python::arg("axis")),(boost::python::arg("positive"))),"add a bounding plane (in fact a sphere with very large radius) bounding the spheres along the direction 'axis' (0,1,2), on the 'positive' or negative side.")
 	.def("setState",&TesselationWrapper::setState,((boost::python::arg("state")=0)),"Make the current state of the simulation the initial (0) or final (1) configuration for the definition of displacement increments, use only state=0 if you just want to get  volmumes and porosity. Exclude bodies using the bitmask :yref:`TesselationWrapper::groupMask`.")
 	.def("loadState",&TesselationWrapper::loadState,(boost::python::arg("inputFile")="state",boost::python::arg("state")=0,boost::python::arg("bz2")=true),"Load a file with positions to define state 0 or 1.")
 	.def("saveState",&TesselationWrapper::saveState,(boost::python::arg("outputFile")="state",boost::python::arg("state")=0,boost::python::arg("bz2")=true),"Save a file with positions, can be later reloaded in order to define state 0 or 1.")
 	.def("volume",&TesselationWrapper::Volume,(boost::python::arg("id")=0),"Returns the volume of Voronoi's cell of a sphere.")
 	.def("defToVtk",&TesselationWrapper::defToVtk,(boost::python::arg("outputFile")="def.vtk"),"Write local deformations in vtk format from states 0 and 1.")
 	.def("defToVtkFromStates",&TesselationWrapper::defToVtkFromStates,(boost::python::arg("input1")="state1",boost::python::arg("input2")="state2",boost::python::arg("outputFile")="def.vtk",boost::python::arg("bz2")=true),"Write local deformations in vtk format from state files (since the file format is very special, consider using defToVtkFromPositions if the input files were not generated by TesselationWrapper).")
 	.def("defToVtkFromPositions",&TesselationWrapper::defToVtkFromPositions,(boost::python::arg("input1")="pos1",boost::python::arg("input2")="pos2",boost::python::arg("outputFile")="def.vtk",boost::python::arg("bz2")=false),"Write local deformations in vtk format from positions files (one sphere per line, with x,y,z,rad separated by spaces).")
 	.def("computeVolumes",&TesselationWrapper::computeVolumes,"compute volumes of all Voronoi's cells.")
	.def("calcVolPoroDef",&TesselationWrapper::calcVolPoroDef,(boost::python::arg("deformation")=false),"Return a table with per-sphere computed quantities. Include deformations on the increment defined by states 0 and 1 if deformation=True (make sure to define states 0 and 1 consistently).")
	.def("computeDeformations",&TesselationWrapper::computeDeformations,"compute per-particle deformation. Get it with :yref:`TesselationWrapper::deformation` (id,i,j).")
	.def("deformation",&TesselationWrapper::deformation,(boost::python::arg("id"),boost::python::arg("i"),boost::python::arg("j")),"Get individual components of the particle deformation tensors")
	.def("deformationTensor",&TesselationWrapper::deformationTensor,(boost::python::arg("id")),"Get particle deformation (tensor)")
	.def("testAlphaShape",&TesselationWrapper::testAlphaShape,(boost::python::arg("alpha")=0),"transitory function, testing AlphaShape feature")
	.def("getAlphaFaces",&TesselationWrapper::getAlphaFaces,(boost::python::arg("alpha")=0),"Get the list of alpha faces for a given alpha. If alpha is not specified or null the minimum alpha resulting in a unique connected domain is used")
	.def("getAlphaCaps",&TesselationWrapper::getAlphaCaps,(boost::python::arg("alpha")=0,boost::python::arg("shrinkedAlpha")=0,boost::python::arg("fixedAlpha")=false),"Get the list of area vectors for the polyhedral caps associated to boundary particles ('extended' alpha-contour). If alpha is not specified or null the minimum alpha resulting in a unique connected domain is used. Taking a smaller 'shrinked' alpha for placing the virtual spheres moves the enveloppe outside the packing, It should be ~(alpha-refRad) typically.")
	.def("applyAlphaForces",&TesselationWrapper::applyAlphaForces,(boost::python::arg("stress"),boost::python::arg("alpha")=0,boost::python::arg("shrinkedAlpha")=0,boost::python::arg("fixedAlpha")=false, boost::python::arg("reset")=true),"set permanent forces based on stress using an alpha shape")
	.def("applyAlphaVel",&TesselationWrapper::applyAlphaVel,(boost::python::arg("velGrad"),boost::python::arg("alpha")=0,boost::python::arg("shrinkedAlpha")=0,boost::python::arg("fixedAlpha")=false),"set velocities based on a velocity gradient tensor using an alpha shape")
	.def("calcAlphaStress",&TesselationWrapper::calcAlphaStress,(boost::python::arg("alpha")=0,boost::python::arg("shrinkedAlpha")=0,boost::python::arg("fixedAlpha")=false),"get the Love-Weber average of the Cauchy stress on the polyhedral caps associated to boundary particles")
	.def("getAlphaGraph",&TesselationWrapper::getAlphaGraph,(boost::python::arg("alpha")=0,boost::python::arg("shrinkedAlpha")=0,boost::python::arg("fixedAlpha")=false),"Get the list of area vectors for the polyhedral caps associated to boundary particles ('extended' alpha-contour). If alpha is not specified or null the minimum alpha resulting in a unique connected domain is used")
	.def("getAlphaVertices",&TesselationWrapper::getAlphaVertices,(boost::python::arg("alpha")=0),"Get the list of 'alpha' bounding spheres for a given alpha. If alpha is not specified or null the minimum alpha resulting in a unique connected domain is used. This function is generating a new alpha shape for each call, not to be used intensively.")
	.def("getVolPoroDef",&TesselationWrapper::calcVolPoroDef,(boost::python::arg("deformation")=false),"|ydeprecated| see new name :yref:`TesselationWrapper.calcVolPoroDef`")
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(TesselationWrapper);
//} // namespace CGT

#ifdef YADE_OPENGL

class GlExtra_AlphaGraph : public GlExtraDrawer {
public:
	bool                                             reset;
	bool                                             refreshDisplay;
	Real                                             alpha;
	Real                                             shrinkedAlpha;
	bool                                             fixedAlpha;
	static GLUquadric*                               gluQuadric;
	static int                                       glCylinderList, oneCylinder;
	vector<Eigen::Transform<Real, 3, Eigen::Affine>> rots;
	vector<Real>                                     lengths;
	vector<Vector3r>                                 pos;
	Real                                             getAlpha() const { return alpha; };
	void                                             setAlpha(Real a)
	{
		reset = true;
		alpha = a;
	};
	Real getShrinkedAlpha() const { return shrinkedAlpha; };
	void setShrinkedAlpha(Real a)
	{
		reset         = true;
		shrinkedAlpha = a;
	};
	bool getFixedAlpha() const { return fixedAlpha; };
	void setFixedAlpha(bool a)
	{
		reset      = true;
		fixedAlpha = a;
	};
	void render() override;
	void refresh() { refreshDisplay = true; };

	// #define BREAK_OPENGL

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(GlExtra_AlphaGraph,GlExtraDrawer,"Display the outer surface defined by alpha contour. Add it to qt.Renderer().extraDrawers. If no instance of TesselationWrapper is provided, the functor will create its own. See :ysrc:`scripts/examples/alphaShapes/GlDrawAlpha.py`.",
		((shared_ptr<TesselationWrapper>,tesselationWrapper,shared_ptr<TesselationWrapper>(),,"Associated instance of TesselationWrapper."))
#ifdef BREAK_OPENGL // for some reason QGLViewer hangs when 'segments' has python binding, reason is unknown. See the #ifndef above.
		((vector<Vector3r>, segments,,,"segments describing the alpha contour"))
#endif
		((Vector3r, color,Vector3r(0.3, 0.3, 0.7),,"color"))
		((Real, lineWidth,3,,"lineWidth in pixels"))
		((Real, radius, 0,,"radius of cylinder representation, if null 1/12th of average diameter will be used"))
		((bool, lighting, true,,"lighting of cylinders"))
		((bool, wire, false,,"display as solid cylinders or lines"))
		,/*ctor*/
		alpha=0; shrinkedAlpha=0; fixedAlpha=false; reset=bool(tesselationWrapper); refreshDisplay=false;
		, /*py*/
		.def("refresh",&GlExtra_AlphaGraph::refresh,"Refresh internals. Particularly usefull for correct display after the TesselationWrapper is modified externally, not needed if 'wire'=True")
		.add_property("alpha",&GlExtra_AlphaGraph::getAlpha,&GlExtra_AlphaGraph::setAlpha,"alpha value")
		.add_property("shrinkedAlpha",&GlExtra_AlphaGraph::getShrinkedAlpha,&GlExtra_AlphaGraph::setShrinkedAlpha,"shrinkedAlpha value")
		.add_property("fixedAlpha",&GlExtra_AlphaGraph::getFixedAlpha,&GlExtra_AlphaGraph::setFixedAlpha,"fixedAlpha option")
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(GlExtra_AlphaGraph);
#endif /*OPENGL*/

} // namespace yade

#endif /* YADE_CGAL */
