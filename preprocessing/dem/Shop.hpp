// 2007 © Václav Šmilauer <eudoxos@arcig.cz>

#pragma once

#include <boost/lambda/lambda.hpp>

#include "lib/base/Logging.hpp"
#include "lib/base/Math.hpp"
#include "core/Body.hpp"

#include <lib/base/AliasNamespaces.hpp>
#include <boost/function.hpp>

namespace yade { // Cannot have #include directive inside.

class Scene;
class Body;
class SimpleViscoelasticBodyParameters;
class ViscElMat;
class FrictMat;
class Interaction;

/*! Miscillaneous utility functions which are believed to be generally useful.
 *
 * All data members are methods are static, no instance of Shop is created. It is not serializable either.
 */

class Shop {
public:
	DECLARE_LOGGER;

	//! create default sphere, along with its bound etc.
	static shared_ptr<Body> sphere(Vector3r center, Real radius, shared_ptr<Material> mat);
	//! create default box with everything needed
	static shared_ptr<Body> box(Vector3r center, Vector3r extents, shared_ptr<Material> mat);
	//! create default tetrahedron
	static shared_ptr<Body> tetra(Vector3r v[4], shared_ptr<Material> mat);

	//! return instance of default FrictMat
	static shared_ptr<FrictMat> defaultGranularMat();

	//! Return vector of pairs (center,radius) loaded from a file with numbers inside
	static vector<boost::tuple<Vector3r, Real, int>>
	loadSpheresFromFile(const string& fname, Vector3r& minXYZ, Vector3r& maxXYZ, Vector3r* cellSize = NULL);

	//! Save spheres in the current simulation into a text file
	static void saveSpheresToFile(string fileName);

	//! Compute the total volume of spheres
	static Real getSpheresVolume(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), int mask = -1);

	//! Compute the total mass of spheres
	static Real getSpheresMass(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), int mask = -1);

	//! Compute porosity; volume must be given for aperiodic simulations
	static Real getPorosity(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), Real volume = -1);

	static Real getPorosityAlt();

	//! Compute porosity by dividing given volume into a grid of voxels;
	static Real getVoxelPorosity(
	        const shared_ptr<Scene>& rb = shared_ptr<Scene>(), int resolution = 500, Vector3r start = Vector3r(0, 0, 0), Vector3r end = Vector3r(0, 0, 0));

	//! Estimate timestep based on P-wave propagation speed
	static Real PWaveTimeStep(const shared_ptr<Scene> rb = shared_ptr<Scene>());

	//! Estimate timestep based on Rayleigh-wave propagation speed
	static Real RayleighWaveTimeStep(const shared_ptr<Scene> rb = shared_ptr<Scene>());

	//! return 2d coordinates of a 3d point within plane defined by rotation axis and inclination of spiral, wrapped to the 0th period
	static boost::tuple<Real, Real, Real>
	spiralProject(const Vector3r& pt, Real dH_dTheta, int axis = 2, Real periodStart = std::numeric_limits<Real>::quiet_NaN(), Real theta0 = 0);

	//! Calculate inscribed circle center of trianlge
	static Vector3r inscribedCircleCenter(const Vector3r& v0, const Vector3r& v1, const Vector3r& v2);

	/// Get viscoelastic parameters kn,cn,ks,cs from analytical solution of
	/// a problem of interaction of pair spheres with mass 1, collision
	/// time tc and restitution coefficients en,es.
	static void getViscoelasticFromSpheresInteraction(Real tc, Real en, Real es, shared_ptr<ViscElMat> b);

	//! Get unbalanced force of the whole simulation
	static Real unbalancedForce(bool useMaxForce = false, Scene* _rb = NULL);
	static Real kineticEnergy(Scene* _rb = NULL, Body::id_t* maxId = NULL);
	//! get total momentum of current simulation
	static Vector3r momentum();
	//! get total angular momentum of current simulation
	static Vector3r angularMomentum(Vector3r origin = Vector3r::Zero());

	static Vector3r totalForceInVolume(Real& avgIsoStiffness, Scene* _rb = NULL);

	//! apply force on contact point on both bodies (reversed on body 2)
	static void applyForceAtContactPoint(
	        const Vector3r& force, const Vector3r& contPt, Body::id_t id1, const Vector3r& pos1, Body::id_t id2, const Vector3r& pos2, Scene* scene);

	//! map scalar variable to 1d colorscale
	static Vector3r scalarOnColorScale(Real x, Real xmin = 0., Real xmax = 1.);

	//! wrap floating number periodically to the given range
	static Real periodicWrap(Real x, Real x0, Real x1, long* period = NULL);

	//! Flip cell shear without affecting interactions, such that abs of shear strain is minimal for each shear component
	static Matrix3r flipCell();

	//! Define the exact average stress in each particle from contour integral ("LW" stands for Love-Weber, since this is what the contour integral gives).
	static void     getStressLWForEachBody(vector<Matrix3r>& bStresses);
	static py::list getStressLWForEachBody();

	//! Compute the dynamic stress tensor for each bodies;
	static py::list getDynamicStress();
	static Matrix3r getTotalDynamicStress(Real volume = 0);

	//! Function to compute overall ("macroscopic") stress.
	static Matrix3r getStress(Real volume = 0);
	static Matrix3r getCapillaryStress(Real volume = 0, bool mindlin = false);
	static Matrix3r stressTensorOfPeriodicCell()
	{
		LOG_WARN("Shop::stressTensorOfPeriodicCelli is DEPRECATED: use getStress instead");
		return Shop::getStress();
	}

	//! Compute the stress tensor in each defined cell, return a stress tensor depth profile
	static py::tuple
	getStressProfile(Real volume, int nCell, Real dz, Real zRef, vector<Real> vPartAverageX, vector<Real> vPartAverageY, vector<Real> vPartAverageZ);
	static py::tuple getStressProfile_contact(Real volume, int nCell, Real dz, Real zRef); //same, only contact contribution

	//! Compute average depth profile of particle velocity (x,y,z) and solid volume fraction. The direction may be specified by direction argument (default is z).
	static py::tuple getDepthProfiles(Real vCell, int nCell, Real dz, Real zRef, bool activateCond = false, Real radiusPy = 0, int direction = 2);
	//! Same, but taking into account point particles
	static py::tuple getDepthProfiles_center(Real vCell, int nCell, Real dz, Real zRef, bool activateCond = false, Real radiusPy = 0);
	// Same as getDepthProfile, but allows to specify the slices on which to average
	static py::tuple getSlicedProfiles(
	        Real         vCell,
	        int          nCell,
	        Real         dP,
	        vector<Real> sliceCenters,
	        vector<Real> sliceWidths,
	        Real         refP,
	        Real         refS,
	        int          dirP         = 2,
	        int          dirS         = 1,
	        bool         activateCond = false,
	        Real         radiusPy     = 0,
	        Real         nSimpson     = 50);

	// Get section area of a sliced sphere
	static Real getSphereSection(Real z, Real R, Real infS, Real supS);

	//! Compute overall ("macroscopic") stress of periodic cell, returning 2 tensors
	//! (contribution of normal and shear forces)
	static py::tuple normalShearStressTensors(bool compressionPositive = false, bool splitNormalTensor = false, Real thresholdForce = NaN);

	//! Function to compute fabric tensor
	static void fabricTensor(
	        Real&            Fmean,
	        Matrix3r&        fabric,
	        Matrix3r&        fabricStrong,
	        Matrix3r&        fabricWeak,
	        Real             cutoff         = 0.0,
	        bool             splitTensor    = false,
	        Real             thresholdForce = NaN,
	        vector<Vector3r> extrema        = vector<Vector3r>());
	static py::tuple fabricTensor(Real cutoff = 0.0, bool splitTensor = false, Real thresholdForce = NaN, vector<Vector3r> extrema = vector<Vector3r>());

	//! Function to set translational and rotational velocities of all bodies to zero
	static void calm(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), int mask = -1);

	//! Get a list of body-ids, which contacts the given body;
	static py::list getBodyIdsContacts(Body::id_t bodyID = -1);

	//! Set material and contact friction to the given value, non-dynamic bodies are not affected
	static void setContactFriction(Real angleRad);

	//! Homothetic change of sizes of spheres and clumps
	static void growParticles(Real multiplier, bool updateMass, bool dynamicOnly);
	//! Change of size of a single sphere or a clump
	// DEPREC, update wrt growParticles()
	static void growParticle(Body::id_t bodyID, Real multiplier, bool updateMass);

	/* \todo implement groupMask */
	static pair<Vector3r, Vector3r> aabbExtrema(Real cutoff = 0.0, bool centers = false);

	//! evaluation of 2D quantities
	static Real getSpheresVolume2D(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), int mask = -1);
	static Real getVoidRatio2D(const shared_ptr<Scene>& rb = shared_ptr<Scene>(), Real zlen = 1);
	//! get stress tensor and tangent operator tensor for FEMxDEM coupling. By Ning Guo
	static py::tuple getStressAndTangent(Real volume = 0, bool symmetry = true);

	//! tests whether p lies in the (bbMin,bbMax) axis-aligned bounding box
	static bool isInBB(Vector3r p, Vector3r bbMin, Vector3r bbMax);
};

} // namespace yade
