/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Anton Gladky, gladky.anton@gmail.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */


/*
 * The class is based on the following article:
 * Molecular-dynamics force models for better control of energy dissipation in numerical simulations of dense granular media
 * Phys. Rev. E 65, 011302 (2001)
 * http://pre.aps.org/abstract/PRE/v65/i1/e011302
 */ 
 
#ifdef PAIR_CLASS

PairStyle(gran/hooke/history/viscel,PairGranHookeHistoryViscEl)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_VISCEL_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_VISCEL_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHookeHistoryViscEl : public PairGranHookeHistory {

 friend class FixWallGranHookeHistorySimple;

 public:

  PairGranHookeHistoryViscEl(class LAMMPS *);
  ~PairGranHookeHistoryViscEl();

  virtual void settings(int, char **);
  virtual void init_granular();

 protected:
  virtual void allocate_properties(int);
  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &,double &);

  //stiffness and damp parameters
  class FixPropertyGlobal *tc1,*e_n1,*e_t1;
  double **k_n,**k_t,**gamma_n,**gamma_t;
  double **tc,**e_n,**e_t;
  
  int damp_massflag;
};

}

#endif
#endif
