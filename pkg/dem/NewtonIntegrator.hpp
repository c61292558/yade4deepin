/*************************************************************************
 Copyright (C) 2008 by Bruno Chareyre		                         *
*  bruno.chareyre@grenoble-inp.fr      					 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include <lib/base/Math.hpp>
#include <core/Callbacks.hpp>
#include <core/Interaction.hpp>
#include <pkg/common/FieldApplier.hpp>
#include <pkg/dem/GlobalStiffnessTimeStepper.hpp>

#ifdef YADE_OPENMP
#include <omp.h>
#endif

namespace yade { // Cannot have #include directive inside.

/*! An engine that can replace the usual series of engines used for integrating the laws of motion.

 */
class State;

class NewtonIntegrator : public FieldApplier {
	inline void     cundallDamp1st(Vector3r& force, const Vector3r& vel);
	inline void     cundallDamp2nd(const Real& dt, const Vector3r& vel, Vector3r& accel);
	inline void     leapfrogTranslate(State*, const Real& dt);                           // leap-frog translate
	inline void     leapfrogSphericalRotate(State*, const Real& dt);                     // leap-frog rotate of spherical body
	inline void     leapfrogAsphericalRotate(State*, const Real& dt, const Vector3r& M); // leap-frog rotate of aspherical body
	inline void     leapfrogAsphericalRotateOmelyan_1998(State*, const Real& dt, const Vector3r& M, int iter);
	inline void     leapfrogAsphericalRotateCarlos_2023(State* state, const Real& dt, const Vector3r& M, int iter);
	Quaternionr     DotQ(const Vector3r& angVel, const Quaternionr& Q);
	inline Vector3r w_dot(const Vector3r w, const Vector3r M, const Vector3r II);

	// compute linear and angular acceleration, respecting State::blockedDOFs
	Vector3r computeAccel(const Vector3r& force, const Real& mass, int blockedDOFs);
	Vector3r computeAngAccel(const Vector3r& torque, const Vector3r& inertia, int blockedDOFs);

	void updateEnergy(const shared_ptr<Body>& b, const State* state, const Vector3r& fluctVel, const Vector3r& f, const Vector3r& m);
#ifdef YADE_OPENMP
	void ensureSync();
	bool syncEnsured;
#endif
	// whether the cell has changed from the previous step
	bool cellChanged;
	int  homoDeform;

	// wether a body has been selected in Qt view
	bool     bodySelected;
	Matrix3r dVelGrad;
	Vector3r dSpin;

	Vector3r computeAccelWithoutGravity(const Vector3r& force, const Real& mass, int blockedDOFs);
	Vector3r addGravity(int blockedDOFs);

public:
	bool densityScaling; // internal for density scaling
	enum class RotAlgorithm { delValle2023 = 1, Omelyan1998 = 2, Fincham1992 = 3 };
	Real updatingDispFactor; //(experimental) Displacement factor used to trigger bound update: the bound is updated only if updatingDispFactor*disp>sweepDist when >0, else all bounds are updated.
	// function to save maximum velocity, for the verlet-distance optimization
	void saveMaximaVelocity(const Body::id_t& id, State* state);
	void saveMaximaDisplacement(const shared_ptr<Body>& b);
	bool get_densityScaling() const;
	void set_densityScaling(bool dsc);

#ifdef YADE_OPENMP
	vector<Real> threadMaxVelocitySq;
#endif
	void action() override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(NewtonIntegrator,GlobalEngine,"Engine integrating newtonian motion equations, see :ref:`sect-motion-integration` for some theoretical background.",
		((Real,damping,0.2,,"damping coefficient for Cundall's non viscous damping (see :ref:`NumericalDamping` and [Chareyre2005]_)"))
		((Vector3r,gravity,Vector3r::Zero(),,"Gravitational acceleration (effectively replaces GravityEngine)."))
		((Real,maxVelocitySq,0,,"stores max. displacement, based on which we trigger collision detection. |yupdate|"))
		((bool,exactAsphericalRot,true,,"Enable more exact body rotation integrator for :yref:`aspherical bodies<Body.aspherical>` *only*, using formulations from [delValle2023]_, [Omelyan1998]_, or [Fincham1992]_ depending on :yref:`rotAlgorithm<NewtonIntegrator.rotAlgorithm>`"))
		((RotAlgorithm,rotAlgorithm,RotAlgorithm::delValle2023,,"Which rotation algorithm to use. Options are: delValle2023, Omelyan1998, Fincham1992."))
		((int,normalizeEvery,5000,,"Normalize the quaternion every normalizeEvery step. Only used in the aspherical formulations from [delValle2023]_, [Omelyan1998]_."))
		((int,niterOmelyan1998,3,,"The number of iterations used to solve the nonlinear system of [Omelyan1998]_ formulation. Provided a small enough timestep, three iterations are enough to make the system converge."))
		((Matrix3r,prevVelGrad,Matrix3r::Zero(),,"Store previous velocity gradient (:yref:`Cell::velGrad`) to track average acceleration in periodic simulations. |yupdate|"))
		#ifdef YADE_BODY_CALLBACK
			((vector<shared_ptr<BodyCallback> >,callbacks,,,"List (std::vector in c++) of :yref:`BodyCallbacks<BodyCallback>` which will be called for each body as it is being processed."))
		#endif
		((Vector3r,prevCellSize,Vector3r(NaN,NaN,NaN),Attr::hidden,":yref:`Cell` size from previous step, used to detect change and find max velocity during periodic simulations"))
		((bool,warnNoForceReset,true,,"Warn when forces were not resetted in this step by :yref:`ForceResetter`; this mostly points to :yref:`ForceResetter` being forgotten incidentally and should be disabled only with a good reason."))
		// energy tracking
		((int,nonviscDampIx,-1,(Attr::hidden|Attr::noSave),"Index of the energy dissipated using the non-viscous damping (:yref:`damping<NewtonIntegrator.damping>`)."))
		((bool,kinSplit,false,,"Whether to separately track translational and rotational kinetic energy."))
                ((bool,dampGravity,true,,"By default, numerical damping applies to ALL forces, even gravity. If this option is set to false, then the gravity forces calculated based on NewtonIntegrator.gravity are excluded from the damping calculation. This option has no effect on gravity forces added by GravityEngine. "))
		((int,kinEnergyIx,-1,(Attr::hidden|Attr::noSave),"Index for kinetic energy in scene->energies."))
		((int,kinEnergyTransIx,-1,(Attr::hidden|Attr::noSave),"Index for translational kinetic energy in scene->energies."))
		((int,kinEnergyRotIx,-1,(Attr::hidden|Attr::noSave),"Index for rotational kinetic energy in scene->energies."))
		((int,mask,-1,,"If mask defined and the bitwise AND between mask and body`s groupMask gives 0, the body will not move/rotate. Velocities and accelerations will be calculated not paying attention to this parameter."))
		,
		/*ctor*/
			timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
			densityScaling=false;
			#ifdef YADE_OPENMP
				threadMaxVelocitySq.resize(omp_get_max_threads()); syncEnsured=false;
			#endif
		,/*py*/
		.add_property("densityScaling",&NewtonIntegrator::get_densityScaling,&NewtonIntegrator::set_densityScaling,"if True, then density scaling [Pfc3dManual30]_ will be applied in order to have a critical timestep equal to :yref:`GlobalStiffnessTimeStepper::targetDt` for all bodies. This option makes the simulation unrealistic from a dynamic point of view, but may speedup quasistatic simulations. In rare situations, it could be useful to not set the scalling factor automatically for each body (which the time-stepper does). In such case revert :yref:`GlobalStiffnessTimeStepper.densityScaling` to False.")
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(NewtonIntegrator);

} // namespace yade
