/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_history_viscel.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryViscEl::PairGranHookeHistoryViscEl(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
    k_n = k_t = gamma_n = gamma_t = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryViscEl::~PairGranHookeHistoryViscEl()
{
    memory->destroy(k_n);
    memory->destroy(k_t);
    memory->destroy(gamma_n);
    memory->destroy(gamma_t);

    // do not destroy coeffFrict, coeffRollFrict, cohEnergyDens
    // since are destroyed in ~PairGranHookeHistory()
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeHistoryViscEl::settings(int narg, char **arg) 
{
    PairGranHookeHistory::settings(narg,arg);
    
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistoryViscEl::init_granular()
{
  int max_type = mpg->max_type();
  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties

  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  
  tc1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("tc","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  e_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("en","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  e_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("et","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  
  
  
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  
  double mpi2 = M_PI*M_PI;
  
  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          
          k_n[i][j] = (mpi2 + pow(log(e_n1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
          k_t[i][j] = 2.0/7.0*(mpi2 + pow(log(e_t1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
          gamma_n[i][j] = -2.0*log(e_n1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
          gamma_t[i][j] = -2.0/7.0*log(e_t1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
          
          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeHistoryViscEl::allocate_properties(int size)
{
    memory->destroy(k_n);
    memory->destroy(k_t);
    memory->destroy(gamma_n);
    memory->destroy(gamma_t);

    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);

    memory->destroy(cohEnergyDens);

    memory->create(k_n,size+1,size+1,"kn");
    memory->create(k_t,size+1,size+1,"kt");
    memory->create(gamma_n,size+1,size+1,"gamman");
    memory->create(gamma_t,size+1,size+1,"gammat");

    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
    memory->create(coeffRollFrict,size+1,size+1,"coeffRollFrict");

    memory->create(cohEnergyDens,size+1,size+1,"cohEnergyDens");
}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */

inline void PairGranHookeHistoryViscEl::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu)
{
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    

          
    kn = meff*k_n[itype][jtype];
    kt = meff*k_t[itype][jtype];

    gamman = meff*gamma_n[itype][jtype];
    gammat = meff*gamma_t[itype][jtype];

    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;

    return;
}
