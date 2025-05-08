/*****************************************************************************
*  2023 © DLH van der Haven, dannyvdhaven@gmail.com, University of Cambridge *
*                                                                            *
*  For details, please see van der Haven et al., "A physically consistent    *
*   Discrete Element Method for arbitrary shapes using Volume-interacting    *
*   Level Sets", Comput. Methods Appl. Mech. Engrg., 414 (116165):1-21       *
*   https://doi.org/10.1016/j.cma.2023.116165                                *
*  This project has been financed by Novo Nordisk A/S (Bagsværd, Denmark).   *
*                                                                            *
*  This program is licensed under GNU GPLv2, see file LICENSE for details.   *
*****************************************************************************/

#pragma once

#include <lib/base/openmp-accu.hpp>
#include <core/Dispatching.hpp>
#include <core/GlobalEngine.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/levelSet/VolumeGeom.hpp>

namespace yade { // Cannot have #include directive inside.

class Law2_VolumeGeom_FrictPhys_Elastic : public LawFunctor {
public:
	OpenMPAccumulator<Real> plasticDissipation;
	bool                    go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	Real                    getPlasticDissipation() const;
	void                    initPlasticDissipation(Real initVal = 0);
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_VolumeGeom_FrictPhys_Elastic,LawFunctor,"Contact law for elasticity, scaling exponentially with the overlap volume, with Mohr-Coulomb plastic failure without cohesion.\nThis law implements a volumetric variant of the classical elastic-plastic law from [CundallStrack1979]_ (see also [Pfc3dManual30]_). The normal force is $F_n=\\min(k_n V_{overlap}^a, 0)$ with $a=1$ (linear) as the default and the convention of positive tensile forces. The shear force is $F_s=k_s u_s$, the plasticity condition defines the maximum value of the shear force: $F_s^{\\max}=F_n\\tan(\\phi)$, with $\\phi$ the friction angle.",
		((Real,volumePower,1.0,,"The exponent $a$ on the overlap volume within the contact law. Setting to 0.5 gives a near-linear relationship of force with respect to penetration distance for spheres."))
		((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only useful if another contact law is used as well)."))
		((bool,traceEnergy,false,,"Define the total energy dissipated in plastic slips at all contacts. This will trace only plastic energy in this law, see O.trackEnergy for a more complete energies tracing."))
		((int,plastDissipIx,-1,(Attr::hidden|Attr::noSave),"Index for plastic dissipation (with O.trackEnergy)"))
		,,
		.def("plasticDissipation",&Law2_VolumeGeom_FrictPhys_Elastic::getPlasticDissipation,"Total energy dissipated in plastic slips at all FrictPhys contacts. Computed only if :yref:`Law2_VolumeGeom_FrictPhys_Elastic::traceEnergy` is true.")
		.def("initPlasticDissipation",&Law2_VolumeGeom_FrictPhys_Elastic::initPlasticDissipation,"Initialize cummulated plastic dissipation to a value (0 by default).")
	);
	// clang-format on
	FUNCTOR2D(VolumeGeom, FrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_VolumeGeom_FrictPhys_Elastic);


class Law2_VolumeGeom_ViscoFrictPhys_Elastic : public Law2_VolumeGeom_FrictPhys_Elastic {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_VolumeGeom_ViscoFrictPhys_Elastic,Law2_VolumeGeom_FrictPhys_Elastic,"Law similar to :yref:`Law2_VolumeGeom_FrictPhys_Elastic` with the addition of shear creep at contacts.",
		((bool,shearCreep,false,," "))
		((Real,viscosity,1,," "))
		((Real,creepStiffness,1,," "))
		,,
	);
	// clang-format on
	FUNCTOR2D(VolumeGeom, ViscoFrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_VolumeGeom_ViscoFrictPhys_Elastic);


} // namespace yade
