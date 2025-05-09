/*************************************************************************
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef YADE_CGAL

#include "MicroMacroAnalyser.hpp"
#include <lib/triangulation/KinematicLocalisationAnalyser.hpp>
#include <lib/triangulation/Tenseur3.h>
#include <lib/triangulation/TriaxialState.h>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/dem/TriaxialCompressionEngine.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((MicroMacroAnalyser));
CREATE_LOGGER(MicroMacroAnalyser);

MicroMacroAnalyser::~MicroMacroAnalyser()
{ /*delete analyser;*/
} //no need, its a shared_ptr now...

void MicroMacroAnalyser::postLoad(MicroMacroAnalyser&)
{
	ofile.open(outputFile.c_str(), std::ios::app);
	if (!boost::filesystem::exists(outputFile.c_str())) ofile << "iteration eps1w eps2w eps3w eps11g eps22g eps33g eps12g eps13g eps23g" << endl;
}

void MicroMacroAnalyser::action()
{
	//cerr << "MicroMacroAnalyser::action() (interval="<< interval <<", iteration="<< scene->iter<<")" << endl;
	if (!triaxialCompressionEngine) {
		vector<shared_ptr<Engine>>::iterator itFirst = scene->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = scene->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "TriaxialCompressionEngine") {
				LOG_DEBUG("stress controller engine found");
				triaxialCompressionEngine = YADE_PTR_CAST<TriaxialCompressionEngine>(*itFirst);
			}
		}
		if (!triaxialCompressionEngine) LOG_ERROR("stress controller engine not found");
	}
	if (triaxialCompressionEngine->strain[0] == 0) return; // no deformation yet
	if (!initialized) {
		setState(1, true, false);
		//Check file here again, to make sure we write to the correct file when filename is modified after the scene is loaded
		ofile.open(outputFile.c_str(), std::ios::app);
		if (!boost::filesystem::exists(outputFile.c_str())) ofile << "iteration eps1w eps2w eps3w eps11g eps22g eps33g eps12g eps13g eps23g" << endl;
		initialized = true;
	} else if (scene->iter % interval == 0) {
		setState(2, true, compIncrt);
		if (compDeformation) {
			analyser->computeParticlesDeformation();
			//for (int i=0; i<analyser->ParticleDeformation.size();i++) cerr<< analyser->ParticleDeformation[i]<<endl;
			std::ostringstream oss;
			oss << "deformation" << incrtNumber++ << ".vtk";
			analyser->DefToFile(oss.str().c_str());
		}
		CGT::Tenseur_sym3 epsg(analyser->grad_u_total);
		ofile << scene->iter << analyser->Delta_epsilon(1, 1) << " " << analyser->Delta_epsilon(2, 2) << " " << analyser->Delta_epsilon(3, 3) << " "
		      << epsg(1, 1) << " " << epsg(2, 2) << " " << epsg(3, 3) << " " << epsg(1, 2) << " " << epsg(1, 3) << " " << epsg(2, 3) << endl;
		analyser->SwitchStates();
	}
	//cerr << "ENDOF MicroMacro::action" << endl;
}

void MicroMacroAnalyser::setState(unsigned int state, bool save_states, bool computeIncrement, mask_t mask)
{
	LOG_INFO("MicroMacroAnalyser::setState");
	CGT::TriaxialState& TS = makeState(state, NULL, mask);
	if (state == 2) {
		analyser->Delta_epsilon(3, 3) = analyser->TS1->eps3 - analyser->TS0->eps3;
		analyser->Delta_epsilon(1, 1) = analyser->TS1->eps1 - analyser->TS0->eps1;
		analyser->Delta_epsilon(2, 2) = analyser->TS1->eps2 - analyser->TS0->eps2;
		if (computeIncrement) {
			analyser->SetForceIncrements();
			analyser->SetDisplacementIncrements();
		}
	}
	if (save_states) {
		std::ostringstream oss;
		//oss<<stateFileName<<"_"<<scene->iter;
		oss << stateFileName << "_" << stateNumber++;
		TS.to_file(oss.str().c_str(), /*use bz2?*/ true);
	}
	LOG_DEBUG("ENDOF MicroMacroAnalyser::setState");
}

//Copy simulation data in the triaxialState structure
CGT::TriaxialState& MicroMacroAnalyser::makeState(unsigned int state, const char* filename, mask_t mask)
{
	//  declaration of ‘scene’ shadows a member of ‘yade::MicroMacroAnalyser’ [-Werror=shadow]
	Scene*                     scene2 = Omega::instance().getScene().get();
	shared_ptr<BodyContainer>& bodies = scene2->bodies;
	CGT::TriaxialState*        ts     = 0;
	if (state == 1) ts = analyser->TS0;
	else if (state == 2)
		ts = analyser->TS1;
	else
		LOG_ERROR("state must be 1 or 2, instead of " << state);
	CGT::TriaxialState& TS = *ts;

	TS.reset();
	auto lengthBodies = bodies->size();
	TS.mean_radius    = 0;
	TS.grains.resize(lengthBodies);
	long               Ng = 0;
	vector<Body::id_t> fictiousVtx;

	for (const auto& bi : *bodies) {
		if (not bi->maskOk(mask)) continue;
		const Body::id_t Idg = bi->getId();
		TS.grains[Idg].id    = Idg;
		TS.maxId             = max(TS.maxId, long(Idg));
		if (not dynamic_cast<Sphere*>(bi->shape.get())) {
			TS.grains[Idg].isSphere = false;
			if (!nonSphereAsFictious or fictiousVtx.size() >= 6) {
				TS.grains[Idg].id = -1; // invalidate so they won't be inserted in triangulation
				continue;
			}
			fictiousVtx.push_back(Idg);
		} else { //then it is a sphere (not a wall)
			++Ng;
			TS.grains[Idg].isSphere = true;
			const Sphere* s         = YADE_CAST<Sphere*>(bi->shape.get());
			//const GranularMat* p = YADE_CAST<GranularMat*> ( (bi)->material.get() );
			const Vector3r& pos = bi->state->pos;
			Real            rad = s->radius;

			TS.grains[Idg].sphere = CGT::Sphere(CGT::Point(pos[0], pos[1], pos[2]), rad);
			//    TS.grains[Idg].translation = trans;
			AngleAxisr aa(bi->state->ori);
			Vector3r   rotVec       = aa.axis() * aa.angle();
			TS.grains[Idg].rotation = CGT::CVector(rotVec[0], rotVec[1], rotVec[2]);
			TS.box.base = CGT::Point(min(TS.box.base.x(), pos.x() - rad), min(TS.box.base.y(), pos.y() - rad), min(TS.box.base.z(), pos.z() - rad));
			TS.box.sommet = CGT::Point(
			        max(TS.box.sommet.x(), pos.x() + rad), max(TS.box.sommet.y(), pos.y() + rad), max(TS.box.sommet.z(), pos.z() + rad));
			TS.mean_radius += TS.grains[Idg].sphere.weight();
		}
	}
	TS.mean_radius /= Ng; //rayon moyen
	LOG_INFO(" loaded : " << Ng << " grains with mean radius = " << TS.mean_radius);
	Real FAR = 1e4;
	if (fictiousVtx.size() < 6) {
		unsigned missing = 6 - fictiousVtx.size();
		TS.grains.resize(lengthBodies + missing);
		for (unsigned fv = lengthBodies; fv < lengthBodies + missing; fv++) {
			fictiousVtx.push_back(fv);
			TS.grains[fv].id       = fv;
			TS.grains[fv].isSphere = false;
		}
	}
	if (fictiousVtx.size() == 6) {
		CGT::Point& Pmin                 = TS.box.base;
		CGT::Point& Pmax                 = TS.box.sommet;
		TS.grains[fictiousVtx[0]].sphere = CGT::Sphere(
		        CGT::Point(0.5 * (Pmin.x() + Pmax.x()), Pmin.y() - FAR * (Pmax.x() - Pmin.x()), 0.5 * (Pmax.z() + Pmin.z())),
		        FAR * (Pmax.x() - Pmin.x()));
		TS.grains[fictiousVtx[1]].sphere = CGT::Sphere(
		        CGT::Point(0.5 * (Pmin.x() + Pmax.x()), Pmax.y() + FAR * (Pmax.x() - Pmin.x()), 0.5 * (Pmax.z() + Pmin.z())),
		        FAR * (Pmax.x() - Pmin.x()));
		TS.grains[fictiousVtx[2]].sphere = CGT::Sphere(
		        CGT::Point(Pmin.x() - FAR * (Pmax.y() - Pmin.y()), 0.5 * (Pmax.y() + Pmin.y()), 0.5 * (Pmax.z() + Pmin.z())),
		        FAR * (Pmax.y() - Pmin.y()));
		TS.grains[fictiousVtx[3]].sphere = CGT::Sphere(
		        CGT::Point(Pmax.x() + FAR * (Pmax.y() - Pmin.y()), 0.5 * (Pmax.y() + Pmin.y()), 0.5 * (Pmax.z() + Pmin.z())),
		        FAR * (Pmax.y() - Pmin.y()));
		TS.grains[fictiousVtx[4]].sphere = CGT::Sphere(
		        CGT::Point(0.5 * (Pmin.x() + Pmax.x()), 0.5 * (Pmax.y() + Pmin.y()), Pmin.z() - FAR * (Pmax.y() - Pmin.y())),
		        FAR * (Pmax.y() - Pmin.y()));
		TS.grains[fictiousVtx[5]].sphere = CGT::Sphere(
		        CGT::Point(0.5 * (Pmin.x() + Pmax.x()), 0.5 * (Pmax.y() + Pmin.y()), Pmax.z() + FAR * (Pmax.y() - Pmin.y())),
		        FAR * (Pmax.y() - Pmin.y()));
	} else
		LOG_INFO(" the number of fictious vertices should be 0 or 6 usually");

	InteractionContainer::iterator ii    = scene2->interactions->begin();
	InteractionContainer::iterator iiEnd = scene2->interactions->end();
	for (; ii != iiEnd; ++ii) {
		if ((*ii)->isReal()) {
			CGT::TriaxialState::Contact* c = new CGT::TriaxialState::Contact;
			TS.contacts.push_back(c);
			CGT::TriaxialState::VectorGrain& grains = TS.grains;
			Body::id_t                       id1    = (*ii)->getId1();
			Body::id_t                       id2    = (*ii)->getId2();

			c->grain1 = &(TS.grains[id1]);
			c->grain2 = &(TS.grains[id2]);
			grains[id1].contacts.push_back(c);
			grains[id2].contacts.push_back(c);
			c->normal = CGT::CVector(
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->normal.x(),
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->normal.y(),
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->normal.z());
			//    c->normal = ( grains[id2].sphere.point()-grains[id1].sphere.point() );
			//    c->normal = c->normal/sqrt ( pow ( c->normal.x(),2 ) +pow ( c->normal.y(),2 ) +pow ( c->normal.z(),2 ) );
			c->position = CGT::CVector(
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->contactPoint.x(),
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->contactPoint.y(),
			        (YADE_CAST<ScGeom*>((*ii)->geom.get()))->contactPoint.z());
			//    c->position = 0.5* ( ( grains[id1].sphere.point()-CGAL::ORIGIN ) +
			//          ( grains[id1].sphere.weight() *c->normal ) +
			//          ( grains[id2].sphere.point()-CGAL::ORIGIN ) -
			//          ( grains[id2].sphere.weight() *c->normal ) );
			c->fn              = YADE_CAST<FrictPhys*>(((*ii)->phys.get()))->normalForce.dot((YADE_CAST<ScGeom*>((*ii)->geom.get()))->normal);
			Vector3r fs        = YADE_CAST<FrictPhys*>((*ii)->phys.get())->shearForce;
			c->fs              = CGT::CVector(fs.x(), fs.y(), fs.z());
			c->old_fn          = c->fn;
			c->old_fs          = c->fs;
			c->frictional_work = 0;
		}
	}
	//Save various parameters if triaxialCompressionEngine is defined
	if (!triaxialCompressionEngine) {
		vector<shared_ptr<Engine>>::iterator itFirst = scene2->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = scene2->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "TriaxialCompressionEngine") {
				LOG_DEBUG("stress controller engine found");
				triaxialCompressionEngine = YADE_PTR_CAST<TriaxialCompressionEngine>(*itFirst);
			}
		}
		if (!triaxialCompressionEngine) LOG_INFO("stress controller engine not found");
	}

	if (triaxialCompressionEngine) {
		TS.wszzh   = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_top][1];
		TS.wsxxd   = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_right][0];
		TS.wsyyfa  = triaxialCompressionEngine->stress[triaxialCompressionEngine->wall_front][2];
		TS.eps3    = triaxialCompressionEngine->strain[2];                      //find_parameter("eps3=", Statefile);
		TS.eps1    = triaxialCompressionEngine->strain[0];                      //find_parameter("eps1=", Statefile);
		TS.eps2    = triaxialCompressionEngine->strain[1];                      //find_parameter("eps2=", Statefile);
		TS.haut    = triaxialCompressionEngine->height;                         //find_parameter("haut=", Statefile);
		TS.larg    = triaxialCompressionEngine->width;                          //find_parameter("larg=", Statefile);
		TS.prof    = triaxialCompressionEngine->depth;                          //find_parameter("prof=", Statefile);
		TS.porom   = 0 /*analyser->computeMacroPorosity() crasher?*/;           //find_parameter("porom=", Statefile);
		TS.ratio_f = triaxialCompressionEngine->ComputeUnbalancedForce(scene2); //find_parameter("ratio_f=", Statefile);
	} else
		TS.wszzh = TS.wsxxd = TS.wsyyfa = TS.eps3 = TS.eps1 = TS.eps2 = TS.haut = TS.larg = TS.prof = TS.porom = TS.ratio_f = 0;
	if (filename != NULL) TS.to_file(filename);
	return TS;
}

// const vector<CGT::Tenseur3>& MicroMacroAnalyser::makeDeformationArray(const char* state_file1, const char* state_file0)
// {
// 	return analyser->computeParticlesDeformation(state_file1, state_file0);
// }

} // namespace yade

#endif /* YADE_CGAL */
