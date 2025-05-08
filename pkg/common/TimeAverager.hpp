// 2024 © Benjamin Dedieu <benjamin.dedieu@proton.me>

#pragma once

#include <core/PartialEngine.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

class TimeAverager : public PartialEngine {
public:
	Vector3r         getPos(Body::id_t id) const;
	Vector3r         getVel(Body::id_t id) const;
	Vector3r         getAngVel(Body::id_t id) const;
	Vector3r         getForce(Body::id_t id) const;
	Vector3r         getTorque(Body::id_t id) const;
	Real             getNbContact(Body::id_t id) const;
	Vector3r         getContactForce(Body::id_t id) const;
	Vector3r         getContactTorque(Body::id_t id) const;
	vector<Vector3r> getContactForceField(Body::id_t id) const;
	void             initialization();

	void action() override;

	// clang-format off
    YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(TimeAverager,PartialEngine,"Average data over time for specific sphere identified by ids. Data are position, velocity, angular velocity, global resultant force and torque, resultant force and torque computed from contacts only and contact force field (see description below). The data must be first initialized with its instantaneous value by running the initialization method. Then averaged values are updated at every time steps with a moving average algorythm, until the initialization method is run again.",
        ((bool,computeContactForceField,false,,"Wether to compute and average contact force field at the surface of the particles. The contact force field is obtained by ditributing the contact forces on a grid at the surface of the sphere. The contact forces are spread on each point of the grid, according to the distance between the contact point and the grid point. The algorythm uses a gaussian kernel to smooth the field. If computeContactForceField is true, grid and sigma parameters must be filled in. This can significantly increase the computation time for dense grid or high number of particles."))
        ((vector<Vector3r>,grid,,,"Grid on which to compute the contact force field. Should be a list of 3D coordinates at the surface of the particle (a simple way to generate a well distributed grid at the surface of a sphere is with the Fibonacci lattice method)."))
        ((Real,sigma,,,"Standard deviation of the Gaussian function, which determines how the contact forces are weighted based on their distance from a contact point. It is usually set at the order of the distance between two points in the grid."))

        ,/*ctor*/
        ,/*py*/
        .def("initialization",&TimeAverager::initialization,"Initialize tAccu to zero and the averaged variables to there instantaneous values. Necessary to execute before any simulation run, otherwise it crashes.")
        .def("getPos",&TimeAverager::getPos,"Get averaged position of particle since last initialization")
        .def("getVel",&TimeAverager::getVel,"Get averaged velocity of particle since last initialization")
        .def("getAngVel",&TimeAverager::getAngVel,"Get averaged angular velocity of particle since last initialization")
        .def("getForce",&TimeAverager::getForce,"Get averaged resultant force of particle since last initialization")
        .def("getTorque",&TimeAverager::getTorque,"Get averaged resultant torque of particle since last initialization")
        .def("getNbContact",&TimeAverager::getNbContact,"Get averaged number of contact points on particle since last initialization")
        .def("getContactForce",&TimeAverager::getContactForce,"Get averaged resultant force computed from contact forces on particle since last initialization")
        .def("getContactTorque",&TimeAverager::getContactTorque,"Get averaged resultant torque computed from contact forces on particle since last initialization")
        .def("getContactForceField",&TimeAverager::getContactForceField,"Get averaged contact force field at the surface of the particle since last initialization")
	);
	// clang-format on

private:
	Real tAccu; // Accumulated time since last initialization

	// Averaged quantities updated at each iteration
	// Output in <map> type  : convenient for C++, but hard to retrieve with boost python
	// so we prefer to use functions getPos, getVel, etc... and make these private attributes.
	boost::unordered_map<Body::id_t, Vector3r>         pos;
	boost::unordered_map<Body::id_t, Vector3r>         vel;
	boost::unordered_map<Body::id_t, Vector3r>         angVel;
	boost::unordered_map<Body::id_t, Vector3r>         force;
	boost::unordered_map<Body::id_t, Vector3r>         torque;
	boost::unordered_map<Body::id_t, Real>             nbContact;
	boost::unordered_map<Body::id_t, Vector3r>         contactForce;
	boost::unordered_map<Body::id_t, Vector3r>         contactTorque;
	boost::unordered_map<Body::id_t, vector<Vector3r>> contactForceField;

	// Util method which updates an averaged quantity with a new instantaneous value
	template <typename T> // template method to work with different types of values (e.g. Vector3, Vector<Vector3>, Real)
	T updateAverage(const T& averagedVal, const T& instantaneousVal, const Real& dt) const;

	// Util method to retrive a value from a map
	template <typename T> // template method to work with different types of values (e.g. Vector3, Vector<Vector3>, Real)
	T getValueFromMap(const boost::unordered_map<Body::id_t, T>& mapObject, Body::id_t id) const;

	// Methods to retrieve instantaneous quantities, i.e. values at the current time step
	Real             getInstantNbContact(const shared_ptr<Body>& b) const;
	Vector3r         getInstantContactForce(const shared_ptr<Body>& b) const;
	Vector3r         getInstantContactTorque(const shared_ptr<Body>& b) const;
	vector<Vector3r> getInstantContactForceField(const shared_ptr<Body>& b) const;
};
REGISTER_SERIALIZABLE(TimeAverager);

} // namespace yade
