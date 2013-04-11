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
    capillarFlag  = false;
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
        else if (strcmp(arg[iarg_],"capillarType") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'capillar_type'");
            iarg_++;
            if(strcmp(arg[iarg_],"weigert") == 0)
                capillarType  = Weigert;
            else if(strcmp(arg[iarg_],"willett") == 0)
                capillarType  = Willett;
            else if(strcmp(arg[iarg_],"herminghaus") == 0)
                capillarType  = Herminghaus;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'weigert', 'willett' or 'herminghaus' after keyword 'capillarType'");
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
  
  if (capillarFlag)  {
    Gamma1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Gamma","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    Theta1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("Theta","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    VB1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("VB","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }
  
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
          
    kn = meff*k_n[itype][jtype];
    kt = meff*k_t[itype][jtype];

    gamman = meff*gamma_n[itype][jtype];
    gammat = meff*gamma_t[itype][jtype];

    xmu=coeffFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    
    return;
}
/* ---------------------------------------------------------------------- */

inline bool PairGranHookeHistoryViscEl::breakContact(int &ip, int &jp, double &rsq, int &touch, int &addflag) {
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
    
    double sCrit = (1+0.5*Theta)*pow(VBCapillar[itype][jtype],1/3.0);
    
    
    
    if (s<sCrit) {
      
      double c0 = 0.96;
      double c1 = 1.1;
      
      Eigen::Vector3f normV = Eigen::Vector3f(x[ip][0] - x[jp][0], x[ip][1] - x[jp][1], x[ip][2] - x[jp][2]);
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
        double sinBeta = pow(Vb/((c0*R*R*R*(1+3*s/R)*(1+c1*sin(Theta)))), 1.0/4.0);
        if ((sinBeta>0.0 and sinBeta<1.0) and (Theta > 0.0 and (Theta <M_PI/2.0))) {
          double beta = asin(sinBeta);
          double r1 = (R*(1-cos(beta)) + s/2.0)/(cos(beta+Theta));
          double r2 = R*sin(beta) + r1*(sin(beta+Theta)-1);
          double Pc = Gamma*(1/r1 + 1/r2);
          fC = 2*M_PI*Gamma*R*sin(beta)*sin(beta+Theta) + M_PI*R*R*Pc*sin(beta)*sin(beta);
        } else {
          touch = 0;
          if (not(sinBeta>0.0 and sinBeta<1.0)){
            error->warning(FLERR,"The Beta is in illegal region!");
          }  else if (not(Theta >= 0.0 and (Theta < M_PI/2.0))) {
            error->warning(FLERR,"The Theta is in illegal region!");
          }
          return true;
        }
      }
      else if (capillarType == Willett) {
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
      else if (capillarType == Herminghaus) {
        double sPl = s/sqrt(Vb/R);
        
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
         
        
        fC = 2.0 * M_PI* R * Gamma * cos(Theta)/(1 + 1.05*sPl + 2.5 *sPl * sPl);         // Herminghaus, equation (7)
      }
      
      
      if (fC != 0.0) {
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
      }
      
      return false;
    } else {
      touch = 0;
      return true;
    }
  }
  
  return true;
};
/* ---------------------------------------------------------------------- */
