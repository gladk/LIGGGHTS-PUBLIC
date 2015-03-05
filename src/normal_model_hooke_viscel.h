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
#include <Eigen/Dense>

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
      gamma(NULL),
      theta(NULL),
      vb(NULL),
      critDist(-1),
      R(-1), vbCur(-1), gammaCur(-1), thetaCur(-1),
      explicitFlag(false),
      capillarFlag(false),
      touchFlag(false)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("explicitly", explicitFlag);
      settings.registerOnOff("capillary", capillarFlag);
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
      
      if (capillarFlag) {
        registry.registerProperty("gamma", &MODEL_PARAMS::createGamma);
        registry.connect("gamma", gamma,"model hooke/viscel");
       
        registry.registerProperty("theta", &MODEL_PARAMS::createTheta);
        registry.connect("theta", theta,"model hooke/viscel");
       
        registry.registerProperty("vb", &MODEL_PARAMS::createVb);
        registry.connect("vb", vb,"model hooke/viscel");
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
      
      if (not(touchFlag) and capillarFlag) {
        touchFlag = true;
        if (critDist < 1) {
          R = 2*cdata.radi*cdata.radj/(cdata.radi + cdata.radj);
          vbCur = vb[itype][jtype];
          gammaCur = gamma[itype][jtype];
          thetaCur = theta[itype][jtype];
          const double Vstar = vbCur/(R*R*R);
          const double Sstar = (1+0.5*thetaCur)*(pow(Vstar,1/3.0) + 0.1*pow(Vstar,2.0/3.0)); // [Willett2000], equation (15), use the full-length e.g 2*Sc
          critDist = Sstar*R;
        }
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

    void noCollision(ContactData& cdata, ForceData& i_forces, ForceData& j_forces){
      if (touchFlag and capillarFlag) {
        cdata.has_force_update = true;
        const Eigen::Vector3d dCur = Eigen::Vector3d(cdata.delta[0],cdata.delta[1],cdata.delta[2]);
        const Eigen::Vector3d dCurN = dCur.normalized();
        const double s = dCur.norm() - cdata.radsum;
        if(critDist  > 0 and s < critDist) {
          /*
          * Capillar model from Willet [Willett2000] (analytical solution), but
          * used also in the work of Herminghaus [Herminghaus2005]
          */
          const double Gamma = gammaCur;
          const double Vb = vbCur;
          const double Theta = thetaCur;
          
          const double sPl = (s/2.0)/sqrt(Vb/R); // [Willett2000], equation (sentence after (11)), s - half-separation, so s*2.0
          const double f_star = cos(Theta)/(1 + 2.1*sPl + 10.0 * pow(sPl, 2.0)); // [Willett2000], equation (12)
          const double Fn = f_star * (2*M_PI*R*Gamma); 
          if(cdata.is_wall) {
            const double area_ratio = cdata.area_ratio;
            const double Fn_ = Fn * cdata.area_ratio;
            i_forces.delta_F[0] = -Fn_ * dCurN[0];
            i_forces.delta_F[1] = -Fn_ * dCurN[1];
            i_forces.delta_F[2] = -Fn_ * dCurN[2];
          } else {
            i_forces.delta_F[0] = -Fn * dCurN[0];
            i_forces.delta_F[1] = -Fn * dCurN[1];
            i_forces.delta_F[2] = -Fn * dCurN[2];
    
            j_forces.delta_F[0] = -i_forces.delta_F[0];
            j_forces.delta_F[1] = -i_forces.delta_F[1];
            j_forces.delta_F[2] = -i_forces.delta_F[2];
          }
        } else {
          touchFlag = false;
        }
      }
    }
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** k_n;
    double ** k_t;
    double ** gamma_n;
    double ** gamma_t;
    
    double ** gamma;
    double ** theta;
    double ** vb;
    
    double ** e_n;
    double ** e_t;
    double ** t_c;
    double ** coeffFrict;
    
    bool explicitFlag;
    bool capillarFlag;
    
    bool touchFlag;
    
    double critDist;
    double R;
    double vbCur;
    double gammaCur;
    double thetaCur;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_VISCEL_H_
#endif
