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

/* ----------------------------------------------------------------------
   Contributing author (2012):
   Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_history_viscel.h"
#include "compute_pair_gran_local.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

#include <Eigen/Core>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryViscEl::PairGranHookeHistoryViscEl(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
    k_n = k_t = gamma_n = gamma_t = GammaCapillar = ThetaCapillar = VBCapillar = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistoryViscEl::~PairGranHookeHistoryViscEl()
{
    memory->destroy(k_n);
    memory->destroy(k_t);
    memory->destroy(gamma_n);
    memory->destroy(gamma_t);
    memory->destroy(GammaCapillar);
    memory->destroy(ThetaCapillar);
    memory->destroy(VBCapillar);
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
  
  Gamma1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Gamma","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  Theta1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Theta","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  VB1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("VB","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  
  double mpi2 = M_PI*M_PI;
  
  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
        /*
         *  Bibtex-entry, 
         * https://github.com/yade/trunk/blob/master/doc/references.bib#L154
         * ============================================
         @Article{ Pournin2001,
         title = "Molecular-dynamics force models for better control of energy dissipation in numerical simulations of dense granular media",
         author = "L. Pournin and Th. M. Liebling and A. Mocellin",
         journal = "Phys. Rev. E",
         volume = "65",
         number = "1",
         pages = "011302",
         numpages = "7",
         year = "2001",
         month = "Dec",
         doi = "10.1103/PhysRevE.65.011302",
         publisher = "American Physical Society"
         *============================================
         * pre.aps.org/abstract/PRE/v65/i1/e011302
         * 
         * Used formula (22)
         */
          
          k_n[i][j] = (mpi2 + pow(log(e_n1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
          k_t[i][j] = 2.0/7.0*(mpi2 + pow(log(e_t1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
          gamma_n[i][j] = -2.0*log(e_n1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
          gamma_t[i][j] = -2.0/7.0*log(e_t1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
          
          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          
          
          /*Capillar
           * 
           * 
           */
           
           GammaCapillar[i][j] = Gamma1->compute_array(i-1,j-1);
           ThetaCapillar[i][j] = Theta1->compute_array(i-1,j-1)*M_PI/180.0;
           VBCapillar[i][j] = VB1->compute_array(i-1,j-1);
           
           
           if ((GammaCapillar[i][j]>=0) and ((ThetaCapillar[i][j]>=0.0) and (ThetaCapillar[i][j]<90.0)) and (VBCapillar[i][j] > 0.0)) {
             cohesionflag = 1;
             error->message(FLERR,"Capillar mode IS activated!");
           } else {
             error->message(FLERR,"Capillar mode is NOT activated!");
           }
             
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
    memory->destroy(GammaCapillar);
    memory->destroy(ThetaCapillar);
    memory->destroy(VBCapillar);

    memory->destroy(coeffFrict);

    memory->create(k_n,size+1,size+1,"kn");
    memory->create(k_t,size+1,size+1,"kt");
    memory->create(gamma_n,size+1,size+1,"gamman");
    memory->create(gamma_t,size+1,size+1,"gammat");
    
    memory->create(GammaCapillar,size+1,size+1,"gammacapillar");
    memory->create(ThetaCapillar,size+1,size+1,"thetacapillar");
    memory->create(VBCapillar,size+1,size+1,"vbcapillar");

    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */

inline void PairGranHookeHistoryViscEl::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu, double &vnnr)
{
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
          
    kn = meff*k_n[itype][jtype];
    kt = meff*k_t[itype][jtype];

    gamman = meff*gamma_n[itype][jtype];
    gammat = meff*gamma_t[itype][jtype];

    xmu=coeffFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    
    /*
    char buffer [50]; 
    int n; 
    double ri = atom->radius[ip];
    double rj = atom->radius[jp];
    //n=sprintf (buffer, "ri = %f,  rj = %f", ri, rj); 
    //error->message(FLERR,buffer);
    if (ri != rj) {
        error->all(FLERR,"Only monodisperse medium can be calculated in capillar mode!");
    }
    */
    
    
    //SCritCapillar = 1.0;
    
    return;
}
/* ---------------------------------------------------------------------- */

inline bool PairGranHookeHistoryViscEl::breakContact(int &ip, int &jp, double &rsq, int &touch, int &addflag) {
  if (touch>0) {
    //r is the distance between the sphere's centeres
    double ri = atom->radius[ip];
    double rj = atom->radius[jp];
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double r = sqrt(rsq);
    //char buffer [50]; int n;n=sprintf (buffer, "ri = %f,  rj = %f", ri, rj); error->message(FLERR,buffer);
    
    
    double c0 = 0.96;
    double c1 = 1.1;
    double R = ri;
    double s = r-ri-rj;
    
    double sCrit = (1+0.5*ThetaCapillar[itype][jtype])*pow(VBCapillar[itype][jtype],1/3.0);
    
    if (s<sCrit) {
      
      double **f = atom->f;
      double **x = atom->x;
      
      Eigen::Vector3f normV = Eigen::Vector3f(x[ip][0] - x[jp][0], x[ip][1] - x[jp][1], x[ip][2] - x[jp][2]);
      normV.normalize();
      
      //char buffer [50]; int n;n=sprintf (buffer, "ri = %f,  rj = %f", s, r); error->message(FLERR,buffer);  
      double beta = asin(pow(VBCapillar[itype][jtype]/((c0*R*R*R*(1+3*s/R)*(1+c1*sin(ThetaCapillar[itype][jtype])))), 1.0/4.0));
      double r1 = (R*(1-cos(beta)) + s/2.0)/(cos(beta+ThetaCapillar[itype][jtype]));
      double r2 = R*sin(beta) + r1*(sin(beta+ThetaCapillar[itype][jtype])-1);
      double Pc = GammaCapillar[itype][jtype]*(1/r1 - 1/r2);
      double fC = 2*M_PI*GammaCapillar[itype][jtype]*R*sin(beta)*sin(beta+ThetaCapillar[itype][jtype]) + M_PI*R*R*Pc*sin(beta)*sin(beta);
      
      Eigen::Vector3f fCV = -fC*normV;
      
      if(computeflag)
      {
        f[ip][0] += fCV(0);
        f[ip][1] += fCV(1);
        f[ip][2] += fCV(2);
      };
      
      if (jp < atom->nlocal && computeflag) {
        f[jp][0] -= fCV(0);
        f[jp][1] -= fCV(1);
        f[jp][2] -= fCV(2);
      };
      
      
      //char buffer [50]; int n;n=sprintf (buffer, "s=%f;  sCrit=%f, step = %d", s, sCrit, update->ntimestep); error->message(FLERR,buffer);
      //n=sprintf (buffer, "f [%f, %f, %f]", fCV(0), fCV(1), fCV(2)); error->message(FLERR,buffer);
      
      if(cpl && addflag) {
        char buffer [50]; int s;s=sprintf (buffer, "i=%d j=%d; %f, %f, %f", ip, jp, fCV(0),fCV(1),fCV(2)); error->message(FLERR,buffer);
        cpl->add_pair(ip,jp,fCV(0),fCV(1),fCV(2),0,0,0,0);
      }
      return false;
    } else {
      touch = 0;
      //char buffer [50]; int n;n=sprintf (buffer, "s=%f;  sCrit=%f!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", s, sCrit); error->message(FLERR,buffer);
      return true;
    }
    //error->message(FLERR,"O0000000000000000000000000000000000000000000000000000000000000!");
  }
  
  return false;
}
/* ---------------------------------------------------------------------- */

inline void PairGranHookeHistoryViscEl::addCohesionForce(int &ip, int &jp,double &r, double &Fn_coh) 
{
    
}

