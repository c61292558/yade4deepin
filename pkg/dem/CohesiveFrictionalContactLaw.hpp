/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>     *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <core/Dispatching.hpp>
#include <core/GlobalEngine.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <boost/tuple/tuple.hpp>

namespace yade { // Cannot have #include directive inside.

// The following code was moved from CohFrictMat.hpp
class CohFrictMat : public FrictMat {
public:
	virtual ~CohFrictMat() {};
	/// Serialization
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CohFrictMat,FrictMat,"Material description extending :yref:`FrictMat` with cohesive properties and rotational stiffnesses. For use e.g. with :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
		((bool,isCohesive,true,,"Whether this body can form possibly cohesive interactions (if true and depending on other parameters such as :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys.setCohesionNow`)."))
		((Real,alphaKr,2.0,,"Dimensionless rolling stiffness."))
		((Real,alphaKtw,2.0,,"Dimensionless twist stiffness."))
		((Real,etaRoll,-1.,,"Dimensionless rolling (aka 'bending') strength. If negative, rolling moment will be elastic."))
		((Real,etaTwist,-1.,,"Dimensionless twisting strength. If negative, twist moment will be elastic."))
		((Real,normalCohesion,-1,,"Tensile strength, homogeneous to a pressure. If negative the normal force is purely elastic."))
		((Real,shearCohesion,-1,,"Shear strength, homogeneous to a pressure. If negative the shear force is purely elastic."))
		((bool,fragile,true,,"does cohesion disappear when contact strength is exceeded?"))
		((bool,momentRotationLaw,false,,"Use bending/twisting moment at contact. The contact may have moments only if both bodies have this flag true. See :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment.always_use_moment_law` for details."))
		,
		createIndex();
		);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(CohFrictMat, FrictMat);
};

REGISTER_SERIALIZABLE(CohFrictMat);

// The following code was moved from CohFrictPhys.hpp
class CohFrictPhys : public RotStiffFrictPhys {
public:
	virtual ~CohFrictPhys() = default;
	void SetBreakingState(bool keepRollingFriction = true)
	{
		cohesionBroken = true;
		normalAdhesion = shearAdhesion = rollingAdhesion = twistingAdhesion = 0;
		if (not keepRollingFriction) {
			maxTwistPl = 0;
			maxTwistPl = 0;
		}
	};

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CohFrictPhys,RotStiffFrictPhys,"An :yref:`interaction physics<IPhys>` that extends :yref:`RotStiffFrictPhys` adding a breakable cohesive nature. Used e.g. by :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`.",
		((bool,cohesionDisablesFriction,false,,"is shear strength the sum of friction and adhesion or only adhesion?"))
		((bool,cohesionBroken,true,Attr::readonly,"is cohesion active? Set to false at the creation of a cohesive contact, and set to true when a fragile contact is broken"))
		((bool,fragile,true,,"do cohesion disappear when contact strength is exceeded?"))
		((Real,maxRollPl,0.0,,"Coefficient of rolling friction (negative means elastic)."))
		((Real,maxTwistPl,0.0,,"Coefficient of twisting friction (negative means elastic)."))
		((Real,normalAdhesion,0,,"tensile strength"))
		((Real,shearAdhesion,0,,"cohesive part of the shear strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`)"))
		((Real,rollingAdhesion,0,,"maximum bending moment when normal force is null (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction` and :yref:`CohFrictPhys::maxRollPl`)"))
		((Real,twistingAdhesion,0,,"maximum twisting moment when normal force is null (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction` and :yref:`CohFrictPhys::maxTwistPl`)"))
		((Real,unp,0,,"plastic normal displacement, only used for tensile behaviour and if :yref:`CohFrictPhys::fragile` =false."))
		((Real,unpMax,-1,,"maximum value of plastic normal displacement (counted positively), after that the interaction breaks even if :yref:`CohFrictPhys::fragile` =false. A negative value (i.e. -1) means no maximum."))
		((bool,momentRotationLaw,false,,"set from :yref:`CohFrictMat::momentRotationLaw` in order to possibly use bending/twisting moment at contacts (if true). See :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::always_use_moment_law` for details."))
		((bool,initCohesion,false,,"Initialize the cohesive behaviour with current state as equilibrium state (same as :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesionNow` but acting on only one interaction)"))
		((Real,creep_viscosity,-1,,"creep viscosity [Pa.s/m]."))
		// internal attributes
		((Vector3r,moment_twist,Vector3r(0,0,0),,"Twist moment"))
		((Vector3r,moment_bending,Vector3r(0,0,0),,"Bending moment"))
		,
		createIndex();
	);
	// clang-format on
	/// Indexable
	REGISTER_CLASS_INDEX(CohFrictPhys, RotStiffFrictPhys);
};

REGISTER_SERIALIZABLE(CohFrictPhys);

class Law2_ScGeom6D_CohFrictPhys_CohesionMoment : public LawFunctor {
public:
	OpenMPAccumulator<Real> plasticDissipation;
	Real                    normElastEnergy();
	Real                    shearElastEnergy();
	Real                    bendingElastEnergy();
	Real                    twistElastEnergy();
	Real                    totalElastEnergy();
	Real                    getPlasticDissipation() const;
	void                    initPlasticDissipation(Real initVal = 0);
	void                    setElasticForces(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I, bool computeMoment, Real& Fn, const Real& dt);
	// This function corrects the trial elastic forces and increment energy dissipation,
	//  it returns true if there is active force/torque, false if the particles are detached (fragile interaction, then the interaction is deleted)
	bool checkPlasticity(ScGeom6D* geom, CohFrictPhys* phys, Real& Fn, bool computeMoment);
	void checkConsistency(const shared_ptr<CohFrictPhys> phys, Body::id_t id1, Body::id_t id2)
	        const; // WARN about possible issues with interaction parameters related to cohesion/friction

	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom6D_CohFrictPhys_CohesionMoment,LawFunctor,"Law for linear traction-compression-bending-twisting, with cohesion+friction and Mohr-Coulomb plasticity surface. This law adds adhesion and moments to :yref:`Law2_ScGeom_FrictPhys_CundallStrack`.\n\nThe normal force is (with the convention of positive tensile forces) $F_n=min(k_n*(u_n-u_n^p), a_n)$, with $a_n$ the normal adhesion and $u_n^p$ the plastic part of normal displacement. The shear force is $F_s=k_s*u_s$, the plasticity condition defines the maximum value of the shear force, by default $F_s^{max}=F_n*tan(\\phi)+a_s$, with $\\phi$ the friction angle and $a_s$ the shear adhesion. If :yref:`CohFrictPhys::cohesionDisablesFriction` is True, friction is ignored as long as adhesion is active, and the maximum shear force is only $F_s^{max}=a_s$.\n\nIf the maximum tensile or maximum shear force is reached and :yref:`CohFrictPhys::fragile` =True (default), the cohesive link is broken, and $a_n, a_s$ are set back to zero. If a tensile force is present, the contact is lost, else the shear strength is $F_s^{max}=F_n*tan(\\phi)$. If :yref:`CohFrictPhys::fragile` =False, the behaviour is perfectly plastic, and the shear strength is kept constant.\n\nIf :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::momentRotationLaw` =True, bending and twisting moments are computed using a linear law with moduli respectively $k_t$ and $k_r$, so that the moments are : $M_b=k_b*\\Theta_b$ and $M_t=k_t*\\Theta_t$, with $\\Theta_{b,t}$ the relative rotations between interacting bodies (details can be found in [Bourrier2013]_). The maximum values of the moments depend on constant terms (:yref:`CohFrictPhys::rollingAdhesion` and :yref:`CohFrictPhys::twistingAdhesion`) and on terms which depend on the normal force $F_n$ through the generalized friction coefficients :yref:`CohFrictPhys::maxRollPl` and :yref:`CohFrictPhys::maxTwistPl`; the instantaneous rolling resistance is thus $a - F_n*\\eta$ if $a$ is adhesion and $\\eta$ is the friction coefficient. \n\nCreep at contact is implemented in this law, as defined in [Hassan2010]_. If activated, there is a viscous behaviour of the shear and twisting components, and the evolution of the elastic parts of shear displacement and relative twist is given by $du_{s,e}/dt=-F_s/\\nu_s$ and $d\\Theta_{t,e}/dt=-M_t/\\nu_t$.\n\nFor turning adhesion on or off during a simulation, see :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesion` (subsequently, it is possible to modify the adhesion values directly, e.g. `i.phys.shearAdhesion=...`)",
		((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
		((bool,always_use_moment_law,false,,"If false, compute moments only for cohesive contacts, broken contacts will have only normal and shear forces. If true, compute bending/twisting moments at all contacts and use the frictional coefficients :yref:`CohFrictMat::etaRoll` and :yref:`CohFrictMat::etaTwist` to define the strength of the broken contacts."))
		((bool,shear_creep,false,,"activate creep on the shear force, using :yref:`CohesiveFrictionalContactLaw::creep_viscosity`."))
		((bool,twist_creep,false,,"activate creep on the twisting moment, using :yref:`CohesiveFrictionalContactLaw::creep_viscosity`."))
		((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at contacts. Note that it will not reflect any energy associated to de-bonding, as it may occur for fragile contacts, nor does it include plastic dissipation in traction."))
		((bool,useIncrementalForm,false,,"use the incremental formulation to compute bending and twisting moments. Creep on the twisting moment is not included in such a case."))
		((int,shearDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for shear dissipation (with O.trackEnergy)"))
		((int,normalDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for normal dissipation (with O.trackEnergy)"))
		((int,bendingDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for bending dissipation (with O.trackEnergy)"))
		((int,twistDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for twist dissipation (with O.trackEnergy)"))
		((bool,consistencyCheck,true,(Attr::hidden|Attr::noSave),"Some parameters check are done when True, then it's turned false."))
		((Real,creep_viscosity,1,,"creep viscosity [Pa.s/m]. probably should be moved to Ip2_CohFrictMat_CohFrictMat_CohFrictPhys."))
		,,
		.def("normElastEnergy",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::normElastEnergy,"Compute normal elastic energy.")
		.def("shearElastEnergy",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::shearElastEnergy,"Compute shear elastic energy.")
		.def("bendingElastEnergy",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::bendingElastEnergy,"Compute bending elastic energy.")
		.def("twistElastEnergy",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::twistElastEnergy,"Compute twist elastic energy.")
		.def("elasticEnergy",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::totalElastEnergy,"Compute total elastic energy.")
		.def("plasticDissipation",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::getPlasticDissipation,"Total energy dissipated in plastic slips at all CohFrictPhys contacts. Computed only if :yref:`Law2_ScGeom_FrictPhys_CundallStrack::traceEnergy` is true.")
		.def("initPlasticDissipation",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
		.def("checkConsistency",&Law2_ScGeom6D_CohFrictPhys_CohesionMoment::checkConsistency,(boost::python::arg("ip"),boost::python::arg("id1"),boost::python::arg("id2")),"Runs consistency checks on an interaction physics and warn if some issues are suspected (e.g. if there is finite cohesion on the shear force and pure elasticity on bending moment). Most inconsistencies should occur when interaction properties are modified in a script; they are not supposed to result from the interaction loop, else please report the bug. This function is called automatically by the law functor on the first instance of a cohesive interaction. To call it in a script: `law.checkConsistency(i.phys,i.id1,i.id2)`.")
	);
	// clang-format on
	FUNCTOR2D(ScGeom6D, CohFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom6D_CohFrictPhys_CohesionMoment);


class CohesiveFrictionalContactLaw : public GlobalEngine {
	shared_ptr<Law2_ScGeom6D_CohFrictPhys_CohesionMoment> functor;

public:
	long iter; /// used for checking if new iteration
	void action() override;

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(CohesiveFrictionalContactLaw,GlobalEngine,"[DEPRECATED] Loop over interactions applying :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment` on all interactions.\n\n.. note::\n  Use :yref:`InteractionLoop` and :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment` instead of this class for performance reasons.",
		((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
		((bool,always_use_moment_law,false,,"If false, compute moments only for cohesive contacts, broken contacts will have only normal and shear forces. If true, compute bending/twisting moments at all contacts and use the frictional coefficients :yref:`CohFrictMat::etaRoll` and :yref:`CohFrictMat::etaTwist` to define the strength of the broken contacts. "))
		((bool,shear_creep,false,,"activate creep on the shear force, using :yref:`CohesiveFrictionalContactLaw::creep_viscosity`."))
		((bool,twist_creep,false,,"activate creep on the twisting moment, using :yref:`CohesiveFrictionalContactLaw::creep_viscosity`."))
		((Real,creep_viscosity,false,,"creep viscosity [Pa.s/m]. probably should be moved to Ip2_CohFrictMat_CohFrictMat_CohFrictPhys..."))
	);
	// clang-format on
};

REGISTER_SERIALIZABLE(CohesiveFrictionalContactLaw);

// The following code was moved from Ip2_CohFrictMat_CohFrictMat_CohFrictPhys.hpp
class Ip2_CohFrictMat_CohFrictMat_CohFrictPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	int  cohesionDefinitionIteration;
	// assign cohesion, return false only if typecasting fails
	bool setCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, CohFrictPhys* contactPhysics = nullptr);
	// the python call will go through dyn_cast since the passed pointer is null
	void pySetCohesion(const shared_ptr<Interaction>& interaction, bool cohesive, bool resetDisp)
	{
		bool assigned = setCohesion(interaction, cohesive, nullptr);
		if (assigned and resetDisp)
			YADE_CAST<CohFrictPhys*>(interaction->phys.get())->unp = YADE_CAST<ScGeom6D*>(interaction->geom.get())->penetrationDepth;
	}

	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Ip2_CohFrictMat_CohFrictMat_CohFrictPhys,IPhysFunctor,
		"Generates cohesive-frictional interactions with moments, used in the contact law :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment`. The normal/shear stiffness and friction definitions are the same as in :yref:`Ip2_FrictMat_FrictMat_FrictPhys`, check the documentation there for details.\n\nAdhesions related to the normal and the shear components are calculated from :yref:`CohFrictMat::normalCohesion` ($C_n$) and :yref:`CohFrictMat::shearCohesion` ($C_s$). For particles of size $R_1$,$R_2$ the adhesion will be $a_i=C_i min(R_1,R_2)^2$, $i=n,s$.\n\nTwist and rolling stiffnesses are proportional to the shear stiffness through dimensionless factors alphaKtw and alphaKr, such that the rotational stiffnesses are defined by $k_s \\alpha_i R_1 R_2$, $i=tw,r. The adhesive contributions to rolling and twisting resistance are by default $C_r=\\frac{C_n}{4} min(R_1,R_2)^3$ and $C_{tw}=\\frac{C_s}{2} min(R_1,R_2)^3$ (inspired by the case of beams).",
		((bool,setCohesionNow,false,,"If true, assign cohesion to all existing contacts in current time-step. The flag is turned false automatically, so that assignment is done in the current timestep only."))
		((bool,setCohesionOnNewContacts,false,,"If true, assign cohesion at all new contacts. If false, only existing contacts can be cohesive (also see :yref:`Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesionNow`), and new contacts are only frictional."))
		((shared_ptr<MatchMaker>,normalCohesion,,,"Instance of :yref:`MatchMaker` determining tensile strength"))
		((shared_ptr<MatchMaker>,shearCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the shear strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`)"))
		((shared_ptr<MatchMaker>,rollingCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the rolling strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`). The default is $\\frac{r}{4}R_t$ with $R_t$ the shear strength (inspired by stress in beams with circular cross-section)."))
		((shared_ptr<MatchMaker>,twistingCohesion,,,"Instance of :yref:`MatchMaker` determining cohesive part of the twisting strength (a frictional term might be added depending on :yref:`CohFrictPhys::cohesionDisablesFriction`). The default is $\\frac{r}{2}R_s$ with $R_s$ the shear strength (inspired by stress in beams with circular cross-section)."))
		((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used."))
		,//ctor
		cohesionDefinitionIteration = -1;
		,//py
		.def("setCohesion",&Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::pySetCohesion,(boost::python::arg("interaction"),boost::python::arg("cohesive"),boost::python::arg("resetDisp")),"Bond or un-bond an interaction with cohesion.\n\n  When ``True``, the resulting state is the same as what is obtained by executing an :yref:`InteractionLoop` with the functor's :yref:`setCohesionNow<Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::setCohesionNow>` or the interaction's :yref:`CohFrictPhys::initCohesion` ``True``. It will use the matchmakers if defined. The only difference is that calling this function explicitly will make the contact cohesive even if not both materials have :yref:`CohFrictMat::isCohesive`=``True``.\n\n When ``False``, the resulting state is the same as after breaking a fragile interaction. If `resetDisp` is ``True``, the current distance is taken as the reference for computing normal displacement and normal force.")
		);
	// clang-format on
	FUNCTOR2D(CohFrictMat, CohFrictMat);
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(Ip2_CohFrictMat_CohFrictMat_CohFrictPhys);

} // namespace yade
