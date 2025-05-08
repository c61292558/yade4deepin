// 2022 © Karol Brzeziński <karol.brze@gmail.com>
#pragma once
#include <core/Scene.hpp>
#include <pkg/common/BoundaryController.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.


class VESupportEngine : public BoundaryController {
private:
	bool needsInit = true;
	void init();

public:
	// displacement accumulated in the dashpots [m]
	vector<Vector3r> uc1, uc2;

	void action() override;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(VESupportEngine,BoundaryController,"Engine that constraints given bodies in place (refPos) with a visco-elastic constrain according to the Burgers model. \n\n.. figure:: fig/rheological-scheme.*\n\t:width: 3cm\n\nBurger's rheological scheme with adopted designations. \n\nThe model of applied constraint can be degenerated to simpler models. Passing negative value of the damping coefficient turns off the corresponding dashpot. A negative value of c2, turns off the whole Kelvin-Voigt branch. By default c1=c2=-1, and model is simplified to an elastic boundary condition. Hence, it can be used as Winkler foundation. \n\nPotential applicatons are presented in [Brzezinski2022]_, and examples section (see :ysrc:`examples/viscoelastic-supports/single-element.py`, and :ysrc:`examples/viscoelastic-supports/discrete-foundation.py`",
			((Real,k1,1,,"Stiffness of spring #1 (the one in Maxwell branch) [N/m]"))
			((Real,k2,1,,"Stiffness of spring #2 (the one in Kelvin-Voigt branch) [N/m]"))
			((Real,c1,-1,,"Damping coeff. of dashpot #1 (the one in Maxwell branch). Negative value turns off the dashpot. [N*s/m]"))
			((Real,c2,-1,,"Damping coeff. of dashpot #2 (the one in Kelvin-Voigt branch). Negative value turns off whole Kelvin-Voigt branch. [N*s/m]"))
			((vector<Body::id_t>,bIds,,,"IDs of bodies that should be attached to supports.")),
			/*ctor*/ needsInit=true;
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(VESupportEngine);
} // namespace yade
