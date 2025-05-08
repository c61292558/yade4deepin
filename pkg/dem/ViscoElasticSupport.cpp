// 2022 © Karol Brzeziński <karol.brze@gmail.com>

#include "ViscoElasticSupport.hpp"
#include <core/Scene.hpp>


namespace yade {

YADE_PLUGIN((VESupportEngine));

/************************ ViscoElasticSupportEngine **********************/
void VESupportEngine::init()
{
	needsInit = false;

	assert(bIds.size() > 0);
	uc1.clear();
	uc2.clear();
	for (size_t i = 0; i < bIds.size(); i++) {
		uc1.push_back(Vector3r::Zero());
		uc2.push_back(Vector3r::Zero());
	}
}

void VESupportEngine::action()
{
	if (needsInit) init();
	assert(bIds.size() == uc1.size() && uc1.size() == uc2.size());

	const Real& dt   = scene->dt;
	long        size = bIds.size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
	for (long counter = 0; counter < size; counter++) {
		const auto id = bIds[counter];
		const auto b  = Body::byId(id, scene);
		if (!b) continue;
		const auto disp = b->state->pos - b->state->refPos;
		const auto f    = k1 * (disp - uc1[counter] - uc2[counter]);
		scene->forces.addForce(id, -1 * f);

		// compute future 'viscous part of displacement'
		const auto fII = f - uc2[counter] * k2;
		if (c1 > 0) { uc1[counter] = uc1[counter] + dt * f / c1; }
		if (c2 > 0) { uc2[counter] = uc2[counter] + dt * fII / c2; }
	}
}
} // namespace yade
