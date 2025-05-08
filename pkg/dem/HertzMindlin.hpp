// 2010 © Chiara Modenese <c.modenese@gmail.com>
//
/*
=== HIGH LEVEL OVERVIEW OF MINDLIN ===

Mindlin is a set of classes to include the Hertz-Mindlin formulation for the contact stiffnesses.
The DMT formulation is also considered (for adhesive particles, rigid and small bodies).

*/

#pragma once

#include <core/Dispatching.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>

#include <lib/base/openmp-accu.hpp>
#include <boost/tuple/tuple.hpp>

namespace yade { // Cannot have #include directive inside.

/******************** MindlinPhys *********************************/
class MindlinPhys : public RotStiffFrictPhys {
public:
	virtual ~MindlinPhys() = default;

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(MindlinPhys,RotStiffFrictPhys,"Representation of an interaction of the Hertz-Mindlin type.",
			((Real,kno,0.0,,"Constant value in the formulation of the normal stiffness"))
			((Real,kso,0.0,,"Constant value in the formulation of the tangential stiffness"))
			((Real,maxBendPl,0.0,,"Coefficient to determine the maximum plastic moment to apply at the contact"))

			((Vector3r,normalViscous,Vector3r::Zero(),,"Normal viscous component"))
			((Vector3r,shearViscous,Vector3r::Zero(),,"Shear viscous component"))
			((Vector3r,shearElastic,Vector3r::Zero(),,"Total elastic shear force"))
			((Vector3r,usElastic,Vector3r::Zero(),,"Total elastic shear displacement (only elastic part)"))
			((Vector3r,usTotal,Vector3r::Zero(),,"Total elastic shear displacement (elastic+plastic part)"))

			//((Vector3r,prevNormal,Vector3r::Zero(),,"Save previous contact normal to compute relative rotation"))
			((Vector3r,momentBend,Vector3r::Zero(),,"Artificial bending moment to provide rolling resistance in order to account for some degree of interlocking between particles"))
			((Vector3r,momentTwist,Vector3r::Zero(),,"Artificial twisting moment (no plastic condition can be applied at the moment)"))
			//((Vector3r,dThetaR,Vector3r::Zero(),,"Incremental rolling vector"))

			((Real,radius,NaN,,"Contact radius (only computed with :yref:`Law2_ScGeom_MindlinPhys_Mindlin::calcEnergy`)"))

			//((Real,gamma,0.0,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
			((Real,adhesionForce,0.0,,"Force of adhesion as predicted by DMT"))
			((bool,isAdhesive,false,,"bool to identify if the contact is adhesive, that is to say if the contact force is attractive"))
			((bool,isSliding,false,,"check if the contact is sliding (useful to calculate the ratio of sliding contacts)"))

			// Contact damping ratios
			((Real,betan,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
			((Real,betas,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))
			((Real,beta,0.0,,"Auxiliary parameter used in the viscous damping model of [Mueller2011]_"))

			// temporary
			((Vector3r,prevU,Vector3r::Zero(),,"Previous local displacement; only used with :yref:`Law2_L3Geom_FrictPhys_HertzMindlin`."))
			((Vector2r,Fs,Vector2r::Zero(),,"Shear force in local axes (computed incrementally)"))
#ifdef PARTIALSAT
			((Real,initD,0,,"initial penetration distance, used for crackaperture estimate"))
			((bool,isBroken,0,,"bool to keep a bond flagged as broken (only useful when displacement criteria is used in partial sat for cracked cell estimates) "))
#endif
			,
			createIndex());
	// clang-format on
	REGISTER_CLASS_INDEX(MindlinPhys, RotStiffFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(MindlinPhys);


/******************** Ip2_FrictMat_FrictMat_MindlinPhys *******/
class Ip2_FrictMat_FrictMat_MindlinPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	FUNCTOR2D(FrictMat, FrictMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(
			Ip2_FrictMat_FrictMat_MindlinPhys,IPhysFunctor, 
				R"""(Calculate physical parameters needed to obtain the normal and shear stiffness values according to the Hertz-Mindlin formulation (no slip solution).\\

There are two available viscous damping models for (1) constant and (2) velocity-dependent coefficient of restitution. In both cases, the viscous forces are calculated as $F_{n,viscous}=c_n \cdot v_n$ ($F_{s,viscous}=c_s \cdot v_s$), where $c_n$ ($c_s$) the normal (shear) viscous damping coefficient and $v_n$ ($v_s$) the normal (shear) component of the relative velocity.\\

(1) Constant coefficient of restitution: The normal (shear) viscous damping coefficient is given by $c_n=2 \cdot \beta_n \cdot \sqrt{m_{bar} \cdot k_n}$ ($c_s=2 \cdot \beta_s \cdot \sqrt{m_{bar} \cdot k_s}$),  where $m_{bar}$ the effective mass, $\beta_n$ ($\beta_s$) normal (shear) viscous damping ratios, and $k_{n}=2 \cdot E^* \cdot \sqrt{R^* \cdot \u_N}$ ($k_{s}=8 \cdot G^* \cdot \sqrt{R \cdot u_N}$) the normal (shear) tangential stiffness values, according to the formulations of Hertz and Mindlin, respectively, and $R^*$, $E^*$, $G^*$ the effective radius, elastic and shear moduli of the interacting particles.

The normal (shear) viscous damping coefficient $c_n$ ($c_s$) can be specified either by providing the normal (shear) viscous damping ratio $\beta_n$ ($\beta_s$), which is then assigned directly to :yref:`MindlinPhys.betan` (:yref:`MindlinPhys.betas`), or by defining the normal (shear) coefficient of restitution $e_n$ ($e_s$) in which case the viscous damping ratios are computed using formula (B6) of [Thornton2013]_, written specifically for the Hertz-Mindlin model (no-slip solution) where the end of contact is considered to take place once the normal force is zero and not once the overlap is zero, thus not allowing attractive elastic forces for non-adhesive contacts, as also discussed in [Schwager2007]_.

(2) Velocity-dependent coefficient of restitution: The viscous damping coefficients are given by $c_n=c_s=A \cdot k_n$, where $A$ a dissipative constant. To calculate this constant, the user has to provide a coefficient of restitution ($e_n$) and an impact velocity ($v_n$) corresponding to this $e_n$, as described in [Mueller2011]_.

The following rules apply:
# It is an error to specify both $e_n$ and $\beta_n$ ($e_s$ and $\beta_s$) or both $v_n$ and $\beta_n$.

# If neither $e_n$ nor $\beta_n$ is given, then :yref:`MindlinPhys.betan` will be zero and no viscous damping will be considered.

# If neither $e_s$ nor $\beta_s$ is given, the value of :yref:`Ip2_FrictMat_FrictMat_MindlinPhys.en` is used for :yref:`Ip2_FrictMat_FrictMat_MindlinPhys.es` and the value of :yref:`MindlinPhys.betan` is used for :yref:`MindlinPhys.betas`, respectively.

The $e_n$, $\beta_n$, $e_s$, $\beta_s$, $v_n$ are :yref:`MatchMaker` objects; they can be constructed from float values to always return constant values.

)""",
			((Real,gamma,0.0,,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
			((Real,eta,0.0,,"Coefficient to determine the plastic bending moment"))
			((Real,krot,0.0,,"Rotational stiffness for moment contact law"))
			((Real,ktwist,0.0,,"Torsional stiffness for moment contact law"))
			((shared_ptr<MatchMaker>,en,,,"Normal coefficient of restitution $e_n$."))
			((shared_ptr<MatchMaker>,es,,,"Shear coefficient of restitution $e_s$."))
			((shared_ptr<MatchMaker>,betan,,,"Normal viscous damping ratio $\\beta_n$."))
			((shared_ptr<MatchMaker>,betas,,,"Shear viscous damping ratio $\\beta_s$."))
			((shared_ptr<MatchMaker>,vn,,,"Impact velocity corresponding to the en value to calculate the dissipative constant $An$ used in the viscous damping model of [Mueller2011]_."))
			((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute the friction angle of an interaction. If ``None``, minimum value is used."))
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_MindlinPhys);


class Law2_ScGeom_MindlinPhys_MindlinDeresiewitz : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*) override;
	FUNCTOR2D(ScGeom, MindlinPhys);
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Law2_ScGeom_MindlinPhys_MindlinDeresiewitz,LawFunctor,
			"Hertz-Mindlin contact law with partial slip solution, as described in [Thornton1991]_.",
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
		);
	// clang-format on
};
REGISTER_SERIALIZABLE(Law2_ScGeom_MindlinPhys_MindlinDeresiewitz);

class Law2_ScGeom_MindlinPhys_HertzWithLinearShear : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*) override;
	FUNCTOR2D(ScGeom, MindlinPhys);
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS(Law2_ScGeom_MindlinPhys_HertzWithLinearShear,LawFunctor,
			"Constitutive law for the Hertz formulation (using :yref:`MindlinPhys.kno`) and linear behavior in shear (using :yref:`MindlinPhys.kso` for stiffness and :yref:`FrictPhys.tangensOfFrictionAngle`). \n\n.. note:: No viscosity or damping. If you need those, look at  :yref:`Law2_ScGeom_MindlinPhys_Mindlin`, which also includes non-linear Mindlin shear.",
				((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
				((int,nonLin,0,,"Shear force nonlinearity (the value determines how many features of the non-linearity are taken in account). 1: ks as in HM 2: shearElastic increment computed as in HM 3. granular ratcheting disabled."))
		);
	// clang-format on
};
REGISTER_SERIALIZABLE(Law2_ScGeom_MindlinPhys_HertzWithLinearShear);


/******************** Law2_ScGeom_MindlinPhys_Mindlin *********/
class Law2_ScGeom_MindlinPhys_Mindlin : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	Real normElastEnergy();
	Real adhesionEnergy();

	Real getfrictionDissipation() const;
	Real getshearEnergy() const;
	Real getnormDampDissip() const;
	Real getshearDampDissip() const;
	Real contactsAdhesive();
	Real ratioSlidingContacts();

	FUNCTOR2D(ScGeom, MindlinPhys);
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(Law2_ScGeom_MindlinPhys_Mindlin,LawFunctor,"Constitutive law for the Hertz-Mindlin formulation. It includes non linear elasticity in the normal direction as predicted by Hertz for two non-conforming elastic contact bodies. In the shear direction, instead, it reseambles the simplified case without slip discussed in Mindlin's paper, where a linear relationship between shear force and tangential displacement is provided. Finally, the Mohr-Coulomb criterion is employed to established the maximum friction force which can be developed at the contact. Moreover, it is also possible to include the effect of linear viscous damping through the definition of the parameters $\\beta_{n}$ and $\\beta_{s}$.",
			((bool,includeAdhesion,false,,"bool to include the adhesion force following the DMT formulation. If true, also the normal elastic energy takes into account the adhesion effect."))
			((bool,calcEnergy,false,,"bool to calculate energy terms (shear potential energy, dissipation of energy due to friction and dissipation of energy due to normal and tangential damping)"))
			((bool,includeMoment,false,,"bool to consider rolling resistance (if :yref:`Ip2_FrictMat_FrictMat_MindlinPhys::eta` is 0.0, no plastic condition is applied.)"))
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))
			((bool,nothing,false,,"dummy attribute for declaring preventGranularRatcheting deprecated"))
			//((bool,LinDamp,true,,"bool to activate linear viscous damping (if false, en and es have to be defined in place of betan and betas)"))

			((OpenMPAccumulator<Real>,frictionDissipation,,Attr::noSave,"Energy dissipation due to sliding"))
			((OpenMPAccumulator<Real>,shearEnergy,,Attr::noSave,"Shear elastic potential energy"))
			((OpenMPAccumulator<Real>,normDampDissip,,Attr::noSave,"Energy dissipated by normal damping"))
			((OpenMPAccumulator<Real>,shearDampDissip,,Attr::noSave,"Energy dissipated by tangential damping"))
			, /* deprec */
			((preventGranularRatcheting, nothing,"this value is no longer used, don't define it."))
			, /* init */
			, /* ctor */
			, /* py */
			.def("contactsAdhesive",&Law2_ScGeom_MindlinPhys_Mindlin::contactsAdhesive,"Compute total number of adhesive contacts.")
			.def("ratioSlidingContacts",&Law2_ScGeom_MindlinPhys_Mindlin::ratioSlidingContacts,"Return the ratio between the number of contacts sliding to the total number at a given time.")
			.def("normElastEnergy",&Law2_ScGeom_MindlinPhys_Mindlin::normElastEnergy,"Compute normal elastic potential energy. It handles the DMT formulation if :yref:`Law2_ScGeom_MindlinPhys_Mindlin::includeAdhesion` is set to true.")
// 			.add_property("preventGranularRatcheting",&Law2_ScGeom_MindlinPhys_Mindlin::deprecAttr,&Law2_ScGeom_MindlinPhys_Mindlin::deprecAttr,"Warn that this is no longer used")
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_MindlinPhys_Mindlin);

// The following code was moved from Ip2_FrictMat_FrictMat_MindlinCapillaryPhys.hpp

class MindlinCapillaryPhys : public MindlinPhys {
public:
	int currentIndexes[4]; // used for faster interpolation (stores previous positions in tables)

	virtual ~MindlinCapillaryPhys() {};

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(MindlinCapillaryPhys,MindlinPhys,"Adds capillary physics to Mindlin's interaction physics.",
				((bool,meniscus,false,Attr::readonly,"True when a meniscus with a non-zero liquid volume (:yref:`vMeniscus<MindlinPhys.vMeniscus>`) has been computed for this interaction"))
				((bool,isBroken,false,,"Might be set to true by the user to make liquid bridge inactive (capillary force is zero)"))
				((Real,capillaryPressure,0.,,"Value of the capillary pressure Uc. Defined as Ugas-Uliquid, obtained from :yref:`corresponding Law2 parameter<Law2_ScGeom_CapillaryPhys_Capillarity.capillaryPressure>`"))
				((Real,vMeniscus,0.,,"Volume of the meniscus"))
				((Real,Delta1,0.,,"Defines the surface area wetted by the meniscus on the smallest grains of radius R1 (R1<R2)"))
				((Real,Delta2,0.,,"Defines the surface area wetted by the meniscus on the biggest grains of radius R2 (R1<R2)"))
				((Vector3r,fCap,Vector3r::Zero(),,"Capillary Force produces by the presence of the meniscus. This is the force acting on particle #2"))
				((short int,fusionNumber,0.,,"Indicates the number of meniscii that overlap with this one"))
				,,createIndex();currentIndexes[0]=currentIndexes[1]=currentIndexes[2]=currentIndexes[3]=0;
				,
				);
	// clang-format on
	REGISTER_CLASS_INDEX(MindlinCapillaryPhys, MindlinPhys);
};
REGISTER_SERIALIZABLE(MindlinCapillaryPhys);


class Ip2_FrictMat_FrictMat_MindlinCapillaryPhys : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;

	FUNCTOR2D(FrictMat, FrictMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ip2_FrictMat_FrictMat_MindlinCapillaryPhys,IPhysFunctor, "RelationShips to use with Law2_ScGeom_CapillaryPhys_Capillarity\n\n In these RelationShips all the interaction attributes are computed. \n\n.. warning::\n\tas in the others :yref:`Ip2 functors<IPhysFunctor>`, most of the attributes are computed only once, when the interaction is new.",
	            ((Real,gamma,0.0,,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
				((Real,eta,0.0,,"Coefficient to determine the plastic bending moment"))
				((Real,krot,0.0,,"Rotational stiffness for moment contact law"))
				((Real,ktwist,0.0,,"Torsional stiffness for moment contact law"))
				((shared_ptr<MatchMaker>,en,,,"Normal coefficient of restitution $e_n$."))
				((shared_ptr<MatchMaker>,es,,,"Shear coefficient of restitution $e_s$."))
				((shared_ptr<MatchMaker>,betan,,,"Normal viscous damping ratio $\\beta_n$."))
				((shared_ptr<MatchMaker>,betas,,,"Shear viscous damping ratio $\\beta_s$."))
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_MindlinCapillaryPhys);

} // namespace yade
