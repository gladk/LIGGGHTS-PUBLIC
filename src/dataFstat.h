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
   Anton Gladky (TU Freiberg)
------------------------------------------------------------------------- */
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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
          ar & _VolWater;
          ar & _DistCurr;
          ar & _DistCrit;
        }

  public:
    Eigen::Vector3d _P1, _P2, _Val;
    int _Id1, _Id2;
    double _VolWater, _DistCurr, _DistCrit;
    DataFstat(Eigen::Vector3d P1, Eigen::Vector3d P2, int Id1, int Id2, Eigen::Vector3d Val, double VolWater, double DistCurr, double DistCrit){
      _P1 = P1;
      _P2 = P2;
      _Id1 = Id1;
      _Id2 = Id2;
      _Val = Val;
      _VolWater = VolWater;
      _DistCurr = DistCurr;
      _DistCrit = DistCrit;
    };
    DataFstat(){
      _P1 = Eigen::Vector3d::Zero();
      _P2 = Eigen::Vector3d::Zero();
      _Id1 = 0.;
      _Id2 = 0.;
      _Val = Eigen::Vector3d::Zero();
      _VolWater = 0;
      _DistCurr = 0;
      _DistCrit = 0;
    };
    void swap_ids_if_needed() {
      if (_Id2 > _Id1) {
        std::swap(_P1, _P2);
        std::swap(_Id1, _Id2);
        _Val*=-1;
      }
    }
};

namespace boost {
  namespace serialization {

    template<class Archive>
      void serialize(Archive & ar, Eigen::Vector3d & g, const unsigned int version)
      {
            ar & g[0];
            ar & g[1];
            ar & g[2];
      }

  }
}
