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
   Contributing author (2012):
   Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
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
#include <Eigen/Core>
#include <vector>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace LAMMPS_NS {

class PairGranHookeHistoryViscEl : public PairGranHookeHistory {

 friend class FixWallGranHookeHistoryViscEl;

 public:

  PairGranHookeHistoryViscEl(class LAMMPS *);
  ~PairGranHookeHistoryViscEl();

  virtual void settings(int, char **);
  virtual void init_granular();

 protected:
  virtual void allocate_properties(int);
  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &, double &, double &);
  
  Eigen::Vector3d breakContact(int &, int &, double &, int &, int &);

  //stiffness and damp parameters
  class FixPropertyGlobal *tc1,*e_n1,*e_t1;
  double **k_n,**k_t,**gamma_n,**gamma_t;
  
  //capillary parameters
  class FixPropertyGlobal *Gamma1, *Theta1, *VB1, *fstat1;
  class FixPropertyGlobal *knSet, *gnSet, *ksSet, *gsSet;
  double **GammaCapillar, **ThetaCapillar, **VBCapillar;
  
  bool capillarFlag;
  bool explicitFlag;
  enum capillar_types_all { Weigert, WillettN, WillettA };
  capillar_types_all capillarType;
  int fstat;
  
  virtual void compute_force(int eflag, int vflag, int addflag);
  
  int damp_massflag;
  long int _ncalls;

  class DataFstat {
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
          {
            ar & _P1;
            ar & _P2;
            ar & _Val;
            ar & _Id1;
            ar & _Id2;
          }

    public:
      Eigen::Vector3d _P1, _P2, _Val;
      int _Id1, _Id2;
      DataFstat(Eigen::Vector3d P1, Eigen::Vector3d P2, int Id1, int Id2, Eigen::Vector3d Val);
      DataFstat() {};
  };

};

}

namespace boost {
  namespace serialization {

    template<class Archive>
      void serialize(Archive & ar, Eigen::Vector3d & g, const unsigned int version)
      {
            ar & g[0];
            ar & g[1];
            ar & g[2];
      }

  } // namespace serialization
} // namespace boost

#endif
#endif
