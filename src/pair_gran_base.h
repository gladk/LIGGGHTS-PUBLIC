/* ----------------------------------------------------------------------
   LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS® is part of CFDEM®project
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
   the producer of the LIGGGHTS® software and the CFDEM®coupling software
   See http://www.cfdem.com/terms-trademark-policy for details.

   LIGGGHTS® is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef PAIR_GRAN_BASE_H_
#define PAIR_GRAN_BASE_H_

#include "contact_interface.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "fix_contact_property_atom.h"
#include "os_specific.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "dataFstat.h"

#include "granular_pair_style.h"

#include <boost/unordered_set.hpp>
#include <utility>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/copy.hpp>

namespace LIGGGHTS {
using namespace ContactModels;

namespace PairStyles {

using namespace LAMMPS_NS;

template<typename ContactModel>
class Granular : private Pointers, public IGranularPairStyle {
  CollisionData * aligned_cdata;
  ForceData * aligned_i_forces;
  ForceData * aligned_j_forces;
  ContactModel cmodel;

  class FixPropertyGlobal *fstat1;
  long int fstat;

  inline void force_update(double * const f, double * const torque,
      const ForceData & forces) {
    for (int coord = 0; coord < 3; coord++) {
      f[coord] += forces.delta_F[coord];
      torque[coord] += forces.delta_torque[coord];
    }
  }

public:
  Granular(class LAMMPS * lmp, PairGran* parent) : Pointers(lmp),
    aligned_cdata(aligned_malloc<CollisionData>(32)),
    aligned_i_forces(aligned_malloc<ForceData>(32)),
    aligned_j_forces(aligned_malloc<ForceData>(32)),
    cmodel(lmp, parent) {
  }

  virtual ~Granular() {
    aligned_free(aligned_cdata);
    aligned_free(aligned_i_forces);
    aligned_free(aligned_j_forces);
  }

  int64_t hashcode()
  { return cmodel.hashcode(); }

  virtual void settings(int nargs, char ** args) {
    Settings settings(lmp);
    cmodel.registerSettings(settings);
    bool success = settings.parseArguments(nargs, args);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== PAIR SETTINGS ====\n");
      settings.print_all(screen);
      fprintf(screen, "==== PAIR SETTINGS ====\n");

      fprintf(logfile, "==== PAIR SETTINGS ====\n");
      settings.print_all(logfile);
      fprintf(logfile, "==== PAIR SETTINGS ====\n");
    }
#endif

    if(!success) {
      error->all(FLERR,settings.error_message.c_str());
    }
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);
    fstat1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("fstat","property/global","scalar",0,0,force->pair_style));
    fstat = fstat1->compute_scalar();

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");
      force->registry.print_all(screen);
      fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");

      fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
      force->registry.print_all(logfile);
      fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
    }
#endif
  }

  virtual void write_restart_settings(FILE * fp)
  {
    int64_t hashcode = ContactModel::STYLE_HASHCODE;
    fwrite(&hashcode, sizeof(int64_t), 1, fp);
  }

  virtual void read_restart_settings(FILE * fp)
  {
    int me = comm->me;
    int64_t hashcode = -1;
    if(me == 0){
      size_t dummy = fread(&hashcode, sizeof(int64_t), 1, fp);
      UNUSED(dummy);
      // sanity check
      if(hashcode != ContactModel::STYLE_HASHCODE)
        error->all(FLERR,"wrong pair style loaded!");
    }
  }

  double stressStrainExponent()
  {
    return cmodel.stressStrainExponent();
  }

  virtual void compute_force(PairGran * pg, int eflag, int vflag, int addflag)
  {
    if (eflag || vflag)
      pg->ev_setup(eflag, vflag);
    else
      pg->evflag = pg->vflag_fdotr = 0;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    const int newton_pair = force->newton_pair;

    int inum = pg->list->inum;
    int * ilist = pg->list->ilist;
    int * numneigh = pg->list->numneigh;

    int ** firstneigh = pg->list->firstneigh;
    int ** firsttouch = pg->listgranhistory ? pg->listgranhistory->firstneigh : NULL;
    double ** firstshear = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

    const int dnum = pg->dnum();
    const bool store_contact_forces = pg->storeContactForces();
    const int freeze_group_bit = pg->freeze_group_bit();

    // clear data, just to be safe
    memset(aligned_cdata, 0, sizeof(CollisionData));
    memset(aligned_i_forces, 0, sizeof(ForceData));
    memset(aligned_j_forces, 0, sizeof(ForceData));
    aligned_cdata->area_ratio = 1.0;

    CollisionData & cdata = *aligned_cdata;
    ForceData & i_forces = *aligned_i_forces;
    ForceData & j_forces = *aligned_j_forces;
    cdata.is_wall = false;
    cdata.computeflag = pg->computeflag();
    cdata.shearupdate = pg->shearupdate();

    cmodel.beginPass(cdata, i_forces, j_forces);
    
    boost::mpi::environment env;
    boost::mpi::communicator world;

    std::vector<DataFstat>  FstatVector;
    
    // loop over neighbors of my atoms

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const double radi = radius[i];
      int * const touch = firsttouch ? firsttouch[i] : NULL;
      double * const allshear = firstshear ? firstshear[i] : NULL;
      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      cdata.i = i;
      cdata.radi = radi;

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double radj = radius[j];
        const double radsum = radi + radj;

        cdata.j = j;
        cdata.delta[0] = delx;
        cdata.delta[1] = dely;
        cdata.delta[2] = delz;
        cdata.rsq = rsq;
        cdata.radj = radj;
        cdata.radsum = radsum;
        cdata.touch = touch ? &touch[jj] : NULL;
        cdata.contact_history = allshear ? &allshear[dnum*jj] : NULL;

        i_forces.reset();
        j_forces.reset();

        if (rsq < radsum * radsum) {
          const double r = sqrt(rsq);
          const double rinv = 1.0 / r;

          // unit normal vector
          const double enx = delx * rinv;
          const double eny = dely * rinv;
          const double enz = delz * rinv;

          // meff = effective mass of pair of particles
          // if I or J part of rigid body, use body mass
          // if I or J is frozen, meff is other particle
          double mi, mj;
          const int itype = type[i];
          const int jtype = type[j];

          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[itype];
            mj = mass[jtype];
          }
          if (pg->fr_pair()) {
            const double * mass_rigid = pg->mr_pair();
            if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
            if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
          }

          double meff = mi * mj / (mi + mj);
          if (mask[i] & freeze_group_bit)
            meff = mj;
          if (mask[j] & freeze_group_bit)
            meff = mi;

          // copy collision data to struct (compiler can figure out a better way to
          // interleave these stores with the double calculations above.
          cdata.itype = itype;
          cdata.jtype = jtype;
          cdata.r = r;
          cdata.rinv = rinv;
          cdata.meff = meff;
          cdata.mi = mi;
          cdata.mj = mj;
          cdata.en[0]   = enx;
          cdata.en[1]   = eny;
          cdata.en[2]   = enz;
          cdata.v_i     = v[i];
          cdata.v_j     = v[j];
          cdata.omega_i = omega[i];
          cdata.omega_j = omega[j];

          cmodel.collision(cdata, i_forces, j_forces);

          // if there is a collision, there will always be a force
          cdata.has_force_update = true;
          
        
          if (not(update->ntimestep % fstat) and (update->ntimestep != 0)) {
            DataFstat FstatTMP(
              Eigen::Vector3d(x[i][0],x[i][1],x[i][2]), 
              Eigen::Vector3d(x[j][0],x[j][1],x[j][2]), 
              atom->tag[i], atom->tag[j], 
              Eigen::Vector3d(i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2]),
              -1,
              -1,
              -1
            );
            FstatTMP.swap_ids_if_needed();
            FstatVector.push_back(FstatTMP);
          }
        } else {
          // apply force update only if selected contact models have requested it
          cdata.has_force_update = false;
          cmodel.noCollision(cdata, i_forces, j_forces);
          if (not(update->ntimestep % fstat) and (update->ntimestep != 0)) {
            DataFstat FstatRet =  cmodel.getDataFstat(cdata, i_forces, j_forces);
            DataFstat FstatTMP(
              Eigen::Vector3d(x[i][0],x[i][1],x[i][2]), 
              Eigen::Vector3d(x[j][0],x[j][1],x[j][2]), 
              atom->tag[i], atom->tag[j], 
              Eigen::Vector3d(i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2]),
              -1,
              -1,
              -1
            );
            FstatTMP._VolWater = FstatRet._VolWater;
            FstatTMP._DistCurr = FstatRet._DistCurr;
            FstatTMP._DistCrit = FstatRet._DistCrit;
            if (FstatTMP._DistCrit) {
              FstatTMP.swap_ids_if_needed();
              FstatVector.push_back(FstatTMP);
            }
          }
        }
        
        if(cdata.has_force_update) {
          if (cdata.computeflag) {
            force_update(f[i], torque[i], i_forces);

            if(newton_pair || j < nlocal) {
              force_update(f[j], torque[j], j_forces);
            }
          }

          if (pg->cpl() && addflag)
            pg->cpl_add_pair(cdata, i_forces);

          if (pg->evflag)
            pg->ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2],cdata.delta[0],cdata.delta[1],cdata.delta[2]);

          if (store_contact_forces)
          {
            double forces_torques_i[6],forces_torques_j[6];

            if(!pg->fix_contact_forces()->has_partner(i,atom->tag[j]))
            {
                vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
                vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
                pg->fix_contact_forces()->add_partner(i,atom->tag[j],forces_torques_i);
            }
            if(!pg->fix_contact_forces()->has_partner(j,atom->tag[i]))
            {
                vectorCopy3D(j_forces.delta_F,&(forces_torques_j[0]));
                vectorCopy3D(j_forces.delta_torque,&(forces_torques_j[3]));
                pg->fix_contact_forces()->add_partner(j,atom->tag[i],forces_torques_j);
            }
          }
        }
      }
    }

    cmodel.endPass(cdata, i_forces, j_forces);
    
    
    if (not(update->ntimestep % fstat) and (update->ntimestep != 0)) {
      if (world.rank() == 0) {
        const auto timestep = update->ntimestep;
        std::vector<std::vector<DataFstat>> allF;
        boost::mpi::gather(world, FstatVector, allF, 0);
        
        
        long int numbForces = 0;
        
        std::vector<DataFstat> commonDataFstat;
        boost::unordered_set<std::pair<int, int> > uniqIdPaitrs;
        
        // Prepare vector of only unique forces
        
        for (int proc = 0; proc < world.size(); ++proc) {
          for (auto FstatTMP : allF[proc] ) {
            auto tmpPair = std::make_pair (FstatTMP._Id1,FstatTMP._Id2);
            auto search = uniqIdPaitrs.find(tmpPair);
            if(search == uniqIdPaitrs.end()) {
              uniqIdPaitrs.insert(tmpPair);
              numbForces++;
              commonDataFstat.push_back(FstatTMP);
            }
          }
        }
        
        
        
        
        
        std::stringstream ss; ss << timestep;
        std::string filename = (std::string("post/fstat_") + ss.str() + std::string(".txt.bz2"));
        
        std::ofstream fileOut(filename.c_str(), std::ios::out | std::ios::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::output> outStream;
        outStream.push(boost::iostreams::bzip2_compressor()); 
        outStream.push(fileOut);
        std::ostream fstatOut(&outStream);
        
        fstatOut << "ITEM: TIMESTEP " << std::endl;
        fstatOut << timestep << std::endl;
        fstatOut << "# "<< std::endl;
        fstatOut << numbForces << std::endl;
        fstatOut << "ITEM: ENTRIES c_fc[1] c_fc[2] c_fc[3] c_fc[4] c_fc[5] c_fc[6] c_fc[7] c_fc[8] c_fc[9] c_fc[10] c_fc[11] c_fc[12] VolWater DistCurr DistCrit" << std::endl;
        for (auto FstatTMP : commonDataFstat ) {
          fstatOut << std::setprecision(15) <<
            FstatTMP._P1(0)  << " " << FstatTMP._P1(1)  << " " << FstatTMP._P1(2)  << " "  << 
            FstatTMP._P2(0)  << " " << FstatTMP._P2(1)  << " " << FstatTMP._P2(2)  << " "  << 
            FstatTMP._Id1  << " " << FstatTMP._Id2  << " 0 " <<
            FstatTMP._Val(0)  << " " << FstatTMP._Val(1)  << " " << FstatTMP._Val(2)  << " "  << 
            FstatTMP._VolWater  << " " << FstatTMP._DistCurr  << " "<< FstatTMP._DistCrit  << std::endl;
        }
      } else {
        boost::mpi::gather(world, FstatVector, 0);
      }
    }
    if(store_contact_forces)
        pg->fix_contact_forces()->do_forward_comm();
  }
};

}

}
#endif /* PAIR_GRAN_BASE_H_ */
