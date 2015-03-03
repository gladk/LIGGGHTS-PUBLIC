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
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HOOKE_VISCEL,hooke/viscel,5)
#else
#ifndef NORMAL_MODEL_HOOKE_VISCEL_H_
#define NORMAL_MODEL_HOOKE_VISCEL_H_
#include "contact_models.h"
#include "math.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include <cmath>

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_VISCEL> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      k_n(NULL),
      k_t(NULL),
      gamma_n(NULL),
      gamma_t(NULL),
      e_n(NULL),
      e_t(NULL),
      t_c(NULL),
      coeffFrict(NULL),
      explicitFlag(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("explicitly", explicitFlag);
    }

    void connectToProperties(PropertyRegistry & registry) {
      
      if (explicitFlag) {
        registry.registerProperty("k_n", &MODEL_PARAMS::createKn);
        registry.connect("k_n", k_n,"model hooke/viscel");
        
        registry.registerProperty("k_t", &MODEL_PARAMS::createKt);
        registry.connect("k_t", k_t,"model hooke/viscel");
  
        
        registry.registerProperty("gamman", &MODEL_PARAMS::createGamman);
        registry.connect("gamman", gamma_n,"model hooke/viscel");
        
        registry.registerProperty("gammat", &MODEL_PARAMS::createGammat);
        registry.connect("gammat", gamma_t,"model hooke/viscel");
      } else {
        registry.registerProperty("e_n", &MODEL_PARAMS::createEn);
        registry.connect("e_n", e_n,"model hooke/viscel");
        
        registry.registerProperty("e_t", &MODEL_PARAMS::createEt);
        registry.connect("e_t", e_t,"model hooke/viscel");
        
        registry.registerProperty("t_c", &MODEL_PARAMS::createTc);
        registry.connect("t_c", t_c,"model hooke/viscel");
      }
      
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hooke/viscel");
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      const double meff = cdata.meff;
      const double xmu  = coeffFrict[cdata.itype][cdata.jtype];
      const double vrel = sqrt(cdata.vtr1*cdata.vtr1 + cdata.vtr2*cdata.vtr2 + cdata.vtr3*cdata.vtr3);

      double kn, kt, gamman, gammat;
      
      if (explicitFlag) {
        kn = k_n[itype][jtype];
        kt = k_t[itype][jtype];
        gamman = gamma_n[itype][jtype];
        gammat = gamma_t[itype][jtype];
      } else {
        kn = (M_PI*M_PI + std::pow(log(e_n[itype][jtype]),2))/(std::pow(t_c[itype][jtype], 2))*meff;
        kt = 2.0/7.0*(M_PI*M_PI + std::pow(log(e_t[itype][jtype]),2))/(std::pow(t_c[itype][jtype], 2))*meff;
        gamman = -2.0 /t_c[itype][jtype] * std::log(e_n[itype][jtype])*meff;
        gammat = -4.0/7.0 /t_c[itype][jtype] * std::log(e_t[itype][jtype])*meff;
      }
      //
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*cdata.vn;    
      const double Fn_contact = kn*(cdata.radsum-cdata.r);
      double Fn               = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if((Fn<0.0))
      {
          Fn = 0.0;
      }
      
      cdata.Fn = Fn;

      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;
      
      // force normalization
      const double Ft_friction = xmu * fabs(Fn);
      const double Ft_damping = gammat*vrel;     
      double Ft = 0.0;

      if (vrel != 0.0) {
        Ft = min(Ft_friction, Ft_damping) / vrel;
      }

      // tangential force due to tangential velocity damping
      const double Ft1 = -Ft*cdata.vtr1;
      const double Ft2 = -Ft*cdata.vtr2;
      const double Ft3 = -Ft*cdata.vtr3;

      // forces & torques
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];
      
      const double tor1 = (eny*Ft3 - enz*Ft2);
      const double tor2 = (enz*Ft1 - enx*Ft3);
      const double tor3 = (enx*Ft2 - eny*Ft1);
      
      // apply normal force
      if(cdata.is_wall) {
        const double area_ratio = cdata.area_ratio;
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0] + Ft1 * area_ratio;
        i_forces.delta_F[1] = Fn_ * cdata.en[1] + Ft2 * area_ratio;
        i_forces.delta_F[2] = Fn_ * cdata.en[2] + Ft3 * area_ratio;
        
        i_forces.delta_torque[0] = -cdata.cri * tor1 * area_ratio;
        i_forces.delta_torque[1] = -cdata.cri * tor2 * area_ratio;
        i_forces.delta_torque[2] = -cdata.cri * tor3 * area_ratio;
      } else {
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0] + Ft1;
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1] + Ft2;
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2] + Ft3;
        i_forces.delta_torque[0] = -cdata.cri * tor1;
        i_forces.delta_torque[1] = -cdata.cri * tor2;
        i_forces.delta_torque[2] = -cdata.cri * tor3;

        j_forces.delta_F[0] = -i_forces.delta_F[0] - Ft1;
        j_forces.delta_F[1] = -i_forces.delta_F[1] - Ft2;
        j_forces.delta_F[2] = -i_forces.delta_F[2] - Ft3;
        j_forces.delta_torque[0] = -cdata.crj * tor1;
        j_forces.delta_torque[1] = -cdata.crj * tor2;
        j_forces.delta_torque[2] = -cdata.crj * tor3;
      }
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** k_n;
    double ** k_t;
    double ** gamma_n;
    double ** gamma_t;
    
    double ** e_n;
    double ** e_t;
    double ** t_c;
    double ** coeffFrict;
    
    bool explicitFlag;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_VISCEL_H_
#endif
