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

#ifdef FIX_CLASS

FixStyle(property/atom/tracer,FixPropertyAtomTracer)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_TRACER_H
#define LMP_FIX_PROPERTY_ATOM_TRACER_H

#include "fix_property_atom.h"

namespace LAMMPS_NS {

enum
{
   MARKER_DIRAC = 0,
   MARKER_HEAVISIDE = 1
};

class FixPropertyAtomTracer : public FixPropertyAtom {

 public:

  FixPropertyAtomTracer(class LAMMPS *, int, char **, bool parse = true);
  ~FixPropertyAtomTracer();

  virtual void init();
  virtual int setmask();
  void end_of_step();
  double compute_scalar();

 protected:

  int iarg_;

  char *tracer_name_;

  int marker_style_;
  int step_;
  int check_every_;
  bool first_mark_;

  // params for region marker
  int iregion_;
  char *idregion_;

  // counter how many particles have been marked
  int nmarked_last_, nmarked_;

}; //end class

}
#endif
#endif
