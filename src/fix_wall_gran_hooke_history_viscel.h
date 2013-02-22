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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (2012):
   Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/hooke/history/viscel,FixWallGranHookeHistoryViscEl)

#else

#ifndef LMP_FIX_WALL_GRAN_HOOKE_HISTORY_VISCEL_H
#define LMP_FIX_WALL_GRAN_HOOKE_HISTORY_VISCEL_H

#include "fix.h"
#include "fix_wall_gran_hooke_history.h"

namespace LAMMPS_NS {

class FixWallGranHookeHistoryViscEl : public FixWallGranHookeHistory {
 public:
  FixWallGranHookeHistoryViscEl(class LAMMPS *, int, char **);
  void init_granular();

 protected:

  virtual void deriveContactModelParams(int, double, double, double &, double &, double &, double &, double &,double &, double &);
  class FixPropertyGlobal *tc1,*e_n1,*e_t1;
  double **k_n,**k_t,**gamma_n,**gamma_t;
  int damp_massflag;
};

}

#endif
#endif
