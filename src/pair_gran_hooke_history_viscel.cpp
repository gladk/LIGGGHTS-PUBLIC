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
#include <iostream>
#include <fstream>
#include "vector_liggghts.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
PairGranHookeHistoryViscEl::DataFstat::DataFstat(Eigen::Vector3d P1, Eigen::Vector3d P2, int Id1, int Id2, Eigen::Vector3d Val){
  _P1 = P1;
  _P2 = P2;
  _Id1 = Id1;
  _Id2 = Id2;
  _Val = Val;
}
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
    capillarFlag  = false;
    explicitFlag  = false;
    capillarType  = Weigert;
    
    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"capilarity") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'capilarity'");
            iarg_++;
            if(strcmp(arg[iarg_],"off") == 0)
                capillarFlag = false;
            else if(strcmp(arg[iarg_],"on") == 0)
                capillarFlag = true;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'capilarity'");
            iarg_++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"explicitly") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Set parameters explicitly: not enough arguments for 'explicit flag'");
            iarg_++;
            if(strcmp(arg[iarg_],"off") == 0)
                explicitFlag = false;
            else if(strcmp(arg[iarg_],"on") == 0)
                explicitFlag = true;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'explicitly'");
            iarg_++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"capillarType") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'capillar_type'");
            iarg_++;
            if(strcmp(arg[iarg_],"weigert") == 0)
                capillarType  = Weigert;
            else if(strcmp(arg[iarg_],"willettA") == 0)
                capillarType  = WillettA;
            else if(strcmp(arg[iarg_],"willettN") == 0)
                capillarType  = WillettN;
            else if(strcmp(arg[iarg_],"rabinovich") == 0)
                capillarType  = Rabinovich;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'weigert', 'willettA', 'willettR',  or 'rabinovich' after keyword 'capillarType'");
            iarg_++;
            hasargs = true;
        } else
        {
            error->all(FLERR,"unknown keyword");
        }
    }       
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
  
  //*fstat*********************************************************

  fstat1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("fstat","property/global","scalar",max_type,max_type,force->pair_style));
  
  //*fstat*********************************************************

  if (capillarFlag)  {
    Gamma1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Gamma","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    Theta1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Theta","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    VB1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("VB","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }
  
  if (explicitFlag)  {
    knSet=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Kn","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    gnSet=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Gn","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    ksSet=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Ks","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    gsSet=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Gs","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }
  
  double mpi2 = M_PI*M_PI;
  fstat = fstat1->compute_scalar();
  
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
          
          
          if (explicitFlag){
            error->message(FLERR,"Explicitly set contact parameters.");
            k_n[i][j] = knSet->compute_array(i-1,j-1);
            k_t[i][j] = ksSet->compute_array(i-1,j-1);
            gamma_n[i][j] = gnSet->compute_array(i-1,j-1);
            gamma_t[i][j] = gsSet->compute_array(i-1,j-1);
          } else {
            if (tc1->compute_array(i-1,j-1) > 0.0) {
              k_n[i][j] = (mpi2 + pow(log(e_n1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
              k_t[i][j] = 2.0/7.0*(mpi2 + pow(log(e_t1->compute_array(i-1,j-1)),2))/(pow(tc1->compute_array(i-1,j-1), 2));
              gamma_n[i][j] = -2.0*log(e_n1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
              gamma_t[i][j] = -2.0/7.0*log(e_t1->compute_array(i-1,j-1))/tc1->compute_array(i-1,j-1);
            } else {
              k_n[i][j] = 0.0;
              k_t[i][j] = 0.0;
              gamma_n[i][j] = 0.0;
              gamma_t[i][j] = 0.0;
            }
          }
            
          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          
          
          /*Capillar
           * 
           * 
           */
           
           if (capillarFlag)  {
             GammaCapillar[i][j] = Gamma1->compute_array(i-1,j-1);
             ThetaCapillar[i][j] = Theta1->compute_array(i-1,j-1)*M_PI/180.0;
             VBCapillar[i][j] = VB1->compute_array(i-1,j-1);
           
           
             if ((GammaCapillar[i][j]>=0) and ((ThetaCapillar[i][j]>=0.0) and (ThetaCapillar[i][j]<90.0)) and (VBCapillar[i][j] > 0.0)) {
               capillarFlag = true;
               error->message(FLERR,"Capillar mode ACTIVATED!");
             } else {
               capillarFlag = false;
               error->message(FLERR,"Capillar mode DISABLED!");
             }
           } else {
             error->message(FLERR,"Capillar mode DISABLED!");
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
          
    if (explicitFlag){
      kn = k_n[itype][jtype];
      kt = k_t[itype][jtype];
  
      gamman = gamma_n[itype][jtype];
      gammat = gamma_t[itype][jtype];
    } else {
      kn = meff*k_n[itype][jtype];
      kt = meff*k_t[itype][jtype];
  
      gamman = meff*gamma_n[itype][jtype];
      gammat = meff*gamma_t[itype][jtype];
    }
    
    xmu=coeffFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    
    return;
}
/* ---------------------------------------------------------------------- */

Eigen::Vector3d PairGranHookeHistoryViscEl::breakContact(int &ip, int &jp, double &rsq, int &touch, int &addflag) {
  if (touch>0 and capillarFlag) {
    //r is the distance between the sphere's centeres
    double ri = atom->radius[ip];
    double rj = atom->radius[jp];
    int itype = atom->type[ip];
    int jtype = atom->type[jp];
    double r = sqrt(rsq);
    
    double **f = atom->f;
    double **x = atom->x;
    
    double R = 2 * ri * rj / (ri + rj);
    double s = (r-ri-rj);
    double Theta = ThetaCapillar[itype][jtype];
    double Vb = VBCapillar[itype][jtype];
    double VbS = Vb/(R*R*R);
    double Gamma = GammaCapillar[itype][jtype];
    
    const double Vstar = VBCapillar[itype][jtype]/(R*R*R);
    double sCrit = (1+0.5*Theta)*(pow(Vstar,1/3.0) + 0.1*pow(Vstar,2.0/3.0))*R;  // [Willett2000], equation (15), use the full-length e.g 2*Sc
    
    if (s<sCrit) {
      Eigen::Vector3d normV = Eigen::Vector3d(x[ip][0] - x[jp][0], x[ip][1] - x[jp][1], x[ip][2] - x[jp][2]);
      normV.normalize();
      double fC = 0.0;
      
      if (capillarType == Weigert) {
       /* Capillar model from Weigert
       * http://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291521-4117%28199910%2916:5%3C238::AID-PPSC238%3E3.0.CO;2-E/abstract
       * 
        ﻿@article {PPSC:PPSC238,
        author = {Weigert, Tom and Ripperger, Siegfried},
        title = {Calculation of the Liquid Bridge Volume and Bulk Saturation from the Half-filling Angle},
        journal = {Particle & Particle Systems Characterization},
        volume = {16},
        number = {5},
        publisher = {WILEY-VCH Verlag GmbH},
        issn = {1521-4117},
        url = {http://dx.doi.org/10.1002/(SICI)1521-4117(199910)16:5<238::AID-PPSC238>3.0.CO;2-E},
        doi = {10.1002/(SICI)1521-4117(199910)16:5<238::AID-PPSC238>3.0.CO;2-E},
        pages = {238--242},
        year = {1999},
        }
        * 
       */ 
          
          double Ca = (1.0 + 6.0*s/(R*2.0));    // [Weigert1999], equation (16)
          double Ct = (1.0 + 1.1*sin(Theta));   // [Weigert1999], equation (17)
          
          double beta = asin(pow(Vb/(0.12*Ca*Ct*pow(2.0*R, 3.0)), 1.0/4.0));                     // [Weigert1999], equation (15), against Vb
          
          double r1 = (2.0*R*(1-cos(beta)) + s)/(2.0*cos(beta+Theta));                           // [Weigert1999], equation (5)
          double r2 = R*sin(beta) + r1*(sin(beta+Theta)-1);                                      // [Weigert1999], equation (6)
          
          double Pk = Gamma*(1/r1 - 1/r2);                                                        /* [Weigert1999], equation (22),
                                                                                                   * see also a sentence over the equation
                                                                                                   * "R1 was taken as positive and R2 was taken as negative"
                                                                                                   */ 
                                                                                                   
          fC = M_PI/4.0*pow((2.0*R),2.0)*pow(sin(beta),2.0)*Pk +
             Gamma*M_PI*2.0*R*sin(beta)*sin(beta+Theta);                                // [Weigert1999], equation (21)
      }
      else if (capillarType == WillettN) {
        double Th1 = Theta;
        double Th2 = Th1*Th1;
        double sPl = s/sqrt(Vb/R)/2.0;
        
        /*
         * Willett, equations in Anhang
        */
        /* Capillar model from Willett
         * http://pubs.acs.org/doi/abs/10.1021/la000657y
         * 
          @article{doi:10.1021/la000657y,
          author = {Willett, Christopher D. and Adams, Michael J. and Johnson, Simon A. and Seville, Jonathan P. K.},
          title = {Capillary Bridges between Two Spherical Bodies},
          journal = {Langmuir},
          volume = {16},
          number = {24},
          pages = {9396-9405},
          year = {2000},
          doi = {10.1021/la000657y},
          
          URL = {http://pubs.acs.org/doi/abs/10.1021/la000657y},
          eprint = {http://pubs.acs.org/doi/pdf/10.1021/la000657y}
          }
         */ 
        double f1 = (-0.44507 + 0.050832*Th1 - 1.1466*Th2) + 
                  (-0.1119 - 0.000411*Th1 - 0.1490*Th2) * log(VbS) +
                  (-0.012101 - 0.0036456*Th1 - 0.01255*Th2) *log(VbS)*log(VbS) +
                  (-0.0005 - 0.0003505*Th1 - 0.00029076*Th2) *log(VbS)*log(VbS)*log(VbS);
        
        double f2 = (1.9222 - 0.57473*Th1 - 1.2918*Th2) +
                  (-0.0668 - 0.1201*Th1 - 0.22574*Th2) * log(VbS) +
                  (-0.0013375 - 0.0068988*Th1 - 0.01137*Th2) *log(VbS)*log(VbS);
                  
                  
        double f3 = (1.268 - 0.01396*Th1 - 0.23566*Th2) +
                  (0.198 + 0.092*Th1 - 0.06418*Th2) * log(VbS) +
                  (0.02232 + 0.02238*Th1 - 0.009853*Th2) *log(VbS)*log(VbS) +
                  (0.0008585 + 0.001318*Th1 - 0.00053*Th2) *log(VbS)*log(VbS)*log(VbS);
        
        double f4 = (-0.010703 + 0.073776*Th1 - 0.34742*Th2) +
                  (0.03345 + 0.04543*Th1 - 0.09056*Th2) * log(VbS) +
                  (0.0018574 + 0.004456*Th1 - 0.006257*Th2) *log(VbS)*log(VbS);
  
        double lnFS = f1 - f2*exp(f3*log(sPl) + f4*log(sPl)*log(sPl));
        double FS = exp(lnFS);
        
        fC = FS * 2.0 * M_PI* R * Gamma;
      }
      else if (capillarType == WillettA) {
        double sPl = (s/2.0)/sqrt(Vb/R);                                            // [Willett2000], equation (sentence after (11)), s - half-separation, so s*2.0
        
        /* Capillar model from Herminghaus (Willett)
         * http://www.tandfonline.com/doi/abs/10.1080/00018730500167855
         * 
          @article{doi:10.1080/00018730500167855,
          author = {Herminghaus * , S.},
          title = {Dynamics of wet granular matter},
          journal = {Advances in Physics},
          volume = {54},
          number = {3},
          pages = {221-261},
          year = {2005},
          doi = {10.1080/00018730500167855},
          
          URL = {http://www.tandfonline.com/doi/abs/10.1080/00018730500167855},
          eprint = {http://www.tandfonline.com/doi/pdf/10.1080/00018730500167855}
          }
         */
        const double f_star = cos(Theta)/(1 + 2.1*sPl + 10.0 * pow(sPl, 2.0));  // [Willett2000], equation (12)
        fC = f_star * (2*M_PI*R*Gamma);                                         // [Willett2000], equation (13), against F
      } else if (capillarType == Rabinovich) {
        /* Capillar model from Rabinovich [Rabinov2005]
         * @article{Rabinov2005,
         *   author = "RABINOVICH Yakov I. and ESAYANUR Madhavan S. and MOUDGIL Brij M.",
         *   institution = "Particle Engineering Research Center, University of Florida, USA; Department of Materials Science and Engineering, University of Florida, USA",
         *   title = "Capillary forces between two spheres with a fixed volume liquid bridge : Theory and experiment",
         *   journal = "Langmuir",
         *   year = "2005",
         *   volume = "21",
         *   number = "24",
         *   pages = "10992--10997",
         *   editor = "American Chemical Society",
         *   issn = "0743-7463",
         *   note = "eng",
         *   keywords = "Atomic force microscopy; Equation; Plane surface; Rupture; Thermodynamics; Methodology; Gas pressure; Calculation; Solid; Interaction energy; Geometry; Plate; Separation; Prediction; Powder; Vapor; Capillary condensation; Theory; Liquid; Force",
         *   url = "http://www.refdoc.fr/Detailnotice?idarticle=7435486"
         *  }
         */
        
        const double H = s;
        const double V = Vb;
        
        double alpha = 0.0;
        double dsp = 0.0;
        if (H!=0.0) {
          alpha = sqrt(H/R*(-1+ sqrt(1 + 2.0*V/(M_PI*R*H*H))));                 // [Rabinov2005], equation (A3)
          
          dsp = H/2.0*(-1.0 + sqrt(1.0 + 2.0*V/(M_PI*R*H*H)));                  // [Rabinov2005], equation (20)
        
          fC = -(2*M_PI*R*Gamma*cos(Theta))/(1+(H/(2*dsp))) - 
                2*M_PI*R*Gamma*sin(alpha)*sin(Theta + alpha);                   // [Rabinov2005], equation (19)
        } else {
          fC = -(2*M_PI*R*Gamma*cos(Theta)) - 
                 2*M_PI*R*Gamma*sin(alpha)*sin(Theta + alpha);                  // [Rabinov2005], equation (19)
        }
        
        fC *=-1;
      }
      
      Eigen::Vector3d fCV = -fC*normV;
        if(computeflag)
        {
            f[ip][0] += fCV(0);
            f[ip][1] += fCV(1);
            f[ip][2] += fCV(2);
        }

        if (jp < atom->nlocal && computeflag) {
          f[jp][0] -= fCV(0);
          f[jp][1] -= fCV(1);
          f[jp][2] -= fCV(2);
        } 
      return fCV;
    } else {
      touch = 0;
      return Eigen::Vector3d::Zero();
    }
  }
  
  return Eigen::Vector3d::Zero();
};
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*~~~~~~~~~~~~~~~~Copy-paste from pair_gran_hooke_history.cpp~~~~~~~~~~~~~*/


void PairGranHookeHistoryViscEl::compute_force(int eflag, int vflag,int addflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu,rmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;
  

  boost::mpi::environment env;
  boost::mpi::communicator world;



  //FSTAT**********************************************
  int iproc;																		//proc number
  double fstattmp,time;																//fstattime, time
  long int timestep;					//timestep
  timestep = update->ntimestep;
  //***************************************************
  
  // loop over neighbors of my atoms

  std::vector<PairGranHookeHistoryViscEl::DataFstat>  FstatVector;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;
      
      Eigen::Vector3d fApply = Eigen::Vector3d::Zero();
      Eigen::Vector3d fCap = Eigen::Vector3d::Zero();

      if (rsq >= radsum*radsum) {
        fCap=breakContact(i, j, rsq, touch[jj], addflag);
        // unset non-touching neighbors
        if (touch[jj] and (fCap==Eigen::Vector3d::Zero())) {
          touch[jj] = 0;
          shear = &allshear[dnum_pairgran*jj];
          shear[0] = 0.0;
          shear[1] = 0.0;
          shear[2] = 0.0;
        }
      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        double deltan=radsum-r;
        cri = radi-0.5*deltan;
        crj = radj-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          itype = type[i];
          jtype = type[j];
          mi = mass[itype];
          mj = mass[jtype];
        }
        if (fix_rigid)
        {
           if(body[i] >= 0) mi = masstotal[body[i]];
           if(body[j] >= 0) mj = masstotal[body[j]];
        }

        meff = mi*mj/(mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr);         

        // normal forces = Hookian contact + normal velocity damping

        damp = gamman*vnnr*rsqinv;  
        ccel = kn*(radsum-r)*rinv - damp;
        
        if (cohesionflag) { 
            addCohesionForce(i,j,r,Fn_coh);
            ccel-=Fn_coh*rinv;
        }

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

        shear = &allshear[dnum_pairgran*jj];

        if (shearupdate && computeflag)
        {
            shear[0] += vtr1*dt;
            shear[1] += vtr2*dt;
            shear[2] += vtr3*dt;

            // rotate shear displacements

            rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
            rsht *= rsqinv;
            shear[0] -= rsht*delx;
            shear[1] -= rsht*dely;
            shear[2] -= rsht*delz;
        }

        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +  shear[2]*shear[2]);

        // tangential forces = shear + tangential velocity damping

        fs1 = - (kt*shear[0]);
        fs2 = - (kt*shear[1]);
        fs3 = - (kt*shear[2]);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        // energy loss from sliding or damping
        if (fs > fn) {
            if (shrmag != 0.0) {
                fs1 *= fn/fs;
                fs2 *= fn/fs;
                fs3 *= fn/fs;
                shear[0] = -fs1/kt;
                shear[1] = -fs2/kt;
                shear[2] = -fs3/kt;
            }
            else fs1 = fs2 = fs3 = 0.0;
        }
        else
        {
            fs1 -= (gammat*vtr1);
            fs2 -= (gammat*vtr2);
            fs3 -= (gammat*vtr3);
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);

        // add rolling friction torque
        vectorZeroize3D(r_torque);
        if(rollingflag)
        {
            vectorSubtract3D(omega[i],omega[j],wr_roll);
            wr_rollmag = vectorMag3D(wr_roll);

            if(wr_rollmag > 0.)
            {
                // calculate torque
                reff=radi*radj/(radi+radj);
                vectorScalarMult3D(wr_roll,rmu*kn*deltan*reff/wr_rollmag,r_torque);

                // remove normal (torsion) part of torque
                double rtorque_dot_delta = r_torque[0]*delx + r_torque[1]*dely + r_torque[2]*delz;
                r_torque_n[0] = delx * rtorque_dot_delta * rsqinv;
                r_torque_n[1] = dely * rtorque_dot_delta * rsqinv;
                r_torque_n[2] = delz * rtorque_dot_delta * rsqinv;
                vectorSubtract3D(r_torque,r_torque_n,r_torque);
            }
        }

        fApply = Eigen::Vector3d(fx, fy, fz);
        if(computeflag)
        {
            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;
            torque[i][0] -= cri*tor1 + r_torque[0];
            torque[i][1] -= cri*tor2 + r_torque[1];
            torque[i][2] -= cri*tor3 + r_torque[2];
        }

        if (j < nlocal && computeflag) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= crj*tor1 - r_torque[0];
          torque[j][1] -= crj*tor2 - r_torque[1];
          torque[j][2] -= crj*tor3 - r_torque[2];
        }

        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
      if (not (timestep % fstat) and 
          (fApply!=Eigen::Vector3d::Zero() or 
           fCap!=Eigen::Vector3d::Zero())) {
        if (fApply==Eigen::Vector3d::Zero()) {
          fApply = fCap;
        }
      
        if (not((j >= nlocal and atom->tag[i] <= atom->tag[j]))) {
        PairGranHookeHistoryViscEl::DataFstat FstatTMP(
            Eigen::Vector3d(x[i][0],x[i][1],x[i][2]), 
            Eigen::Vector3d(x[j][0],x[j][1],x[j][2]), 
            atom->tag[i], atom->tag[j], 
            fApply);        
        FstatVector.push_back(FstatTMP);
        }
      }
    }
  }
    
  if ((not (timestep % fstat)) and 
      (world.rank() == 0) and
      (timestep != 0)) {
  
      std::vector< std::vector <PairGranHookeHistoryViscEl::DataFstat > > allF;
      boost::mpi::gather(world, FstatVector, allF, 0);
      
      std::ofstream fstatOut;
      long int numbForces = 0;
      for (int proc = 0; proc < world.size(); ++proc) {
        numbForces += allF[proc].size();
      }
  
      std::string filename;
      std::ostringstream oss;
      oss << "post/fstat_" << timestep <<".txt";
      filename += oss.str();
      fstatOut.open(filename.c_str(), std::ios::out);

      fstatOut << "ITEM: TIMESTEP " << std::endl;
      fstatOut << timestep << std::endl;
      fstatOut << "# "<< std::endl;
      fstatOut << numbForces << std::endl;
      fstatOut << "ITEM: ENTRIES c_fc[1] c_fc[2] c_fc[3] c_fc[4] c_fc[5] c_fc[6] c_fc[7] c_fc[8] c_fc[9] c_fc[10] c_fc[11] c_fc[12] " << std::endl;
      for (int proc = 0; proc < world.size(); ++proc) {
          for (long int i=0; i<allF[proc].size(); i++){
          PairGranHookeHistoryViscEl::DataFstat FstatTMP = allF[proc][i];
            fstatOut << 
              FstatTMP._P1(0)  << " " << FstatTMP._P1(1)  << " " << FstatTMP._P1(2)  << " "  << 
              FstatTMP._P2(0)  << " " << FstatTMP._P2(1)  << " " << FstatTMP._P2(2)  << " "  << 
              FstatTMP._Id1  << " " << FstatTMP._Id2  << " 0 " <<
              FstatTMP._Val(0)  << " " << FstatTMP._Val(1)  << " " << FstatTMP._Val(2)  << " "  << std::endl;
          }
        }
        fstatOut.close();
      } else if (not (timestep % fstat) and (timestep != 0)) {
        boost::mpi::gather(world, FstatVector, 0);
      }
}
