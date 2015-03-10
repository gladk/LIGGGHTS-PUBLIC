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
#include <vector>
#include <Eigen/Dense>


namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_VISCEL> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup) : Pointers(lmp),
      e_n(NULL),
      e_t(NULL),
      t_c(NULL),
      coeffFrict(NULL),
      gamma(NULL),
      theta(NULL),
      vb(NULL),
      capillaryModelLoad(NULL),
      explicitFlag(false),
      capillarFlag(false)
    {
      history_offset = hsetup->add_history_value("firstTouch", "1");
      hsetup->add_history_value("touchFlag", "1");
      hsetup->add_history_value("firstTouchCap", "1");
      hsetup->add_history_value("critDist", "1");
      hsetup->add_history_value("R", "1");
      hsetup->add_history_value("vbCur", "1");
      hsetup->add_history_value("thetaCur", "1");
      hsetup->add_history_value("kn", "1");
      hsetup->add_history_value("kt", "1");
      hsetup->add_history_value("gamman", "1");
      hsetup->add_history_value("gammat", "1");
    }
    
    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("explicitly", explicitFlag);
      settings.registerOnOff("capillary", capillarFlag);
    }
    
    static double Willett_analytic_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      /*
      * Capillar model from Willet [Willett2000] (analytical solution), but
      * used also in the work of Herminghaus [Herminghaus2005]
      */
      
      const double sPl = (s/2.0)/sqrt(Vb/R); // [Willett2000], equation (sentence after (11)), s - half-separation, so s*2.0
      const double f_star = cos(Theta)/(1 + 2.1*sPl + 10.0 * pow(sPl, 2.0)); // [Willett2000], equation (12)
      const double fC = f_star * (2*M_PI*R*Gamma); // [Willett2000], equation (13), against F
      return fC;
    }
    
    static double Willett_numeric_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      const double VbS = Vb/(R*R*R);
      const double Th1 = Theta;
      const double Th2 = Theta;
      
      /*
      * [Willett2000], equations in Attachment
      */
      const double f1 = (-0.44507 + 0.050832*Th1 - 1.1466*Th2) +
        (-0.1119 - 0.000411*Th1 - 0.1490*Th2) * log(VbS) +
        (-0.012101 - 0.0036456*Th1 - 0.01255*Th2) *log(VbS)*log(VbS) +
        (-0.0005 - 0.0003505*Th1 - 0.00029076*Th2) *log(VbS)*log(VbS)*log(VbS);
      const double f2 = (1.9222 - 0.57473*Th1 - 1.2918*Th2) +
        (-0.0668 - 0.1201*Th1 - 0.22574*Th2) * log(VbS) +
        (-0.0013375 - 0.0068988*Th1 - 0.01137*Th2) *log(VbS)*log(VbS);
      const double f3 = (1.268 - 0.01396*Th1 - 0.23566*Th2) +
        (0.198 + 0.092*Th1 - 0.06418*Th2) * log(VbS) +
        (0.02232 + 0.02238*Th1 - 0.009853*Th2) *log(VbS)*log(VbS) +
        (0.0008585 + 0.001318*Th1 - 0.00053*Th2) *log(VbS)*log(VbS)*log(VbS);
      const double f4 = (-0.010703 + 0.073776*Th1 - 0.34742*Th2) +
        (0.03345 + 0.04543*Th1 - 0.09056*Th2) * log(VbS) +
        (0.0018574 + 0.004456*Th1 - 0.006257*Th2) *log(VbS)*log(VbS);
      const double sPl = (s/2.0)/sqrt(Vb/R);
      const double lnFS = f1 - f2*exp(f3*log(sPl) + f4*log(sPl)*log(sPl));
      const double FS = exp(lnFS);
      const double fC = FS * 2.0 * M_PI* R * Gamma;
      return fC;
    }
    
    static double Weigert_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      /*
      * Capillar model from [Weigert1999]
      */
      
      const double a = s;
      const double Ca = (1.0 + 6.0*a/(R*2.0)); // [Weigert1999], equation (16)
      const double Ct = (1.0 + 1.1*sin(Theta)); // [Weigert1999], equation (17)
      
      const double beta = asin(pow(Vb/(0.12*Ca*Ct*pow(2.0*R, 3.0)), 1.0/4.0)); // [Weigert1999], equation (15), against Vb
      const double r1 = (2.0*R*(1-cos(beta)) + a)/(2.0*cos(beta+Theta)); // [Weigert1999], equation (5)
      const double r2 = R*sin(beta) + r1*(sin(beta+Theta)-1); // [Weigert1999], equation (6)
      const double Pk = Gamma*(1/r1 - 1/r2); // [Weigert1999], equation (22),
      // see also a sentence over the equation
      // "R1 was taken as positive and R2 was taken as negative"
      // fC = M_PI*2.0*R*phys.gamma/(1+tan(0.5*beta)); // [Weigert1999], equation (23), [Fisher]
      const double fC = M_PI/4.0*pow((2.0*R),2.0)*pow(sin(beta),2.0)*Pk + Gamma*M_PI*2.0*R*sin(beta)*sin(beta+Theta); // [Weigert1999], equation (21)
      return fC; 
    }
    
    
    static double Rabinovich_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      /*
      * Capillar model from Rabinovich [Rabinov2005]
      *
      * This formulation from Rabinovich has been later verified and corrected
      * by Lambert [Lambert2008]. So we can calculate both formulations
      *
      */
      const double H = s;
      const double V = Vb;
      double fC = 0.0;
      double dsp = 0.0;
      if (H!=0.0) {
        dsp = H/2.0*(-1.0 + sqrt(1.0 + 2.0*V/(M_PI*R*H*H))); // [Rabinov2005], equation (20)
        fC = -(2*M_PI*R*Gamma*cos(Theta))/(1+(H/(2*dsp))); // [Lambert2008], equation (65), taken from [Rabinov2005]
      const double alpha = sqrt(H/R*(-1+ sqrt(1 + 2.0*V/(M_PI*R*H*H)))); // [Rabinov2005], equation (A3)
        fC -= 2*M_PI*R*Gamma*sin(alpha)*sin(Theta + alpha); // [Rabinov2005], equation (19)
      } else {
        fC = -(2*M_PI*R*Gamma*cos(Theta));
        const double alpha = 0.0;
        fC -= 2*M_PI*R*Gamma*sin(alpha)*sin(Theta + alpha); // [Rabinov2005], equation (19)
      }
      fC *=-1;
      return fC;
    }
    
    static double Lambert_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      /*
      * Capillar model from Rabinovich [Rabinov2005]
      *
      * This formulation from Rabinovich has been later verified and corrected
      * by Lambert [Lambert2008]. So we can calculate both formulations
      *
      */
      const double H = s;
      const double V = Vb;
      double fC = 0.0;
      double dsp = 0.0;
      if (H!=0.0) {
        dsp = H/2.0*(-1.0 + sqrt(1.0 + 2.0*V/(M_PI*R*H*H))); // [Rabinov2005], equation (20)
        fC = -(2*M_PI*R*Gamma*cos(Theta))/(1+(H/(2*dsp))); // [Lambert2008], equation (65), taken from [Rabinov2005]
      } else {
        fC = -(2*M_PI*R*Gamma*cos(Theta));
      }
      fC *=-1;
      return fC;
    }
    
    static double Soulie_f(const double & R, const double & Vb, const double & Gamma, const double & Theta, const double & s) {
      /*
      * Capillar model from Soulie [Soulie2006]
      *
      * !!! In this implementation the radiis of particles are taken equal
      * to get the symmetric forces.
      *
      * Please, use this model only for testing purposes.
      *
      */
      
      const double D = s;
      const double V = Vb;
      const double a = -1.1*pow((V/(R*R*R)), -0.53);
      const double b = (-0.148*log(V/(R*R*R)) - 0.96)*Theta*Theta -0.0082*log(V/(R*R*R)) + 0.48;
      const double c = 0.0018*log(V/(R*R*R)) + 0.078;
      const double fC = M_PI*Gamma*sqrt(R*R)*(c+exp(a*D/R+b));
      return fC;  
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
        
        registry.registerProperty("capillaryModel", &MODEL_PARAMS::createCapillaryModel);
        registry.connect("capillaryModel", capillaryModelLoad,"model hooke/viscel");
        
        capModels.push_back(Willett_analytic_f);
        capModels.push_back(Willett_numeric_f);
        capModels.push_back(Weigert_f);
        capModels.push_back(Rabinovich_f);
        capModels.push_back(Lambert_f);
        capModels.push_back(Soulie_f);
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
      
      double * const history = &cdata.contact_history[history_offset];
      double * const firstTouch = &history[0];
      double * const touchFlag = &history[1];
      double * const firstTouchCap = &history[2];
      double * const critDist = &history[3];
      double * const R = &history[4];
      double * const vbCur = &history[5];
      double * const thetaCur = &history[6];
      double * const gammaCur = &history[7];
      double * const kn = &history[8];
      double * const kt = &history[9];
      double * const gamman = &history[10];
      double * const gammat = &history[11];
      
      //std::cout<<"1a  "<<firstTouch<<"  "<<*firstTouch<<"; Time: "<<update->ntimestep<<"; i: "<<cdata.i<<"; j:"<<cdata.j<<std::endl;
      
      if (not(static_cast<int>(*firstTouch))) {
        if (explicitFlag) {
          cdata.kn = k_n[itype][jtype];
          cdata.kt = k_t[itype][jtype];
          cdata.gamman = gamma_n[itype][jtype];
          cdata.gammat = gamma_t[itype][jtype];
        } else {
          cdata.kn = (M_PI*M_PI + std::pow(log(e_n[itype][jtype]),2))/(std::pow(t_c[itype][jtype], 2))*meff;
          cdata.kt = 2.0/7.0*(M_PI*M_PI + std::pow(log(e_t[itype][jtype]),2))/(std::pow(t_c[itype][jtype], 2))*meff;
          cdata.gamman = -2.0 /t_c[itype][jtype] * std::log(e_n[itype][jtype])*meff;
          cdata.gammat = -4.0/7.0 /t_c[itype][jtype] * std::log(e_t[itype][jtype])*meff;
        }
        
        *kn = cdata.kn;
        *kt = cdata.kt;
        *gamman = cdata.gamman;
        *gammat = cdata.gammat;
        
        *firstTouch = 1;
        *touchFlag = 1;
        *firstTouchCap = 1;
      }
      
      cdata.kn = *kn;
      cdata.kt = *kt;
      cdata.gamman = *gamman;
      cdata.gammat = *gammat;
      
      // convert Kn and Kt from pressure units to force/distance^2
      cdata.kn /= force->nktv2p;
      cdata.kt /= force->nktv2p;

      const double Fn_damping = -cdata.gamman*cdata.vn;    
      const double Fn_contact = cdata.kn*(cdata.radsum-cdata.r);
      double Fn               = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if((Fn<0.0))
      {
          Fn = 0.0;
      }
      
      cdata.Fn = Fn;
      
      // force normalization
      const double Ft_friction = xmu * fabs(Fn);
      const double Ft_damping = cdata.gammat*vrel;     
      double Ft = 0.0;

      if (vrel != 0.0) {
        Ft = min(Ft_friction, Ft_damping) / vrel;
      }
      
      if (*firstTouchCap and capillarFlag) {
        *R = 2*cdata.radi*cdata.radj/(cdata.radi + cdata.radj);
        *vbCur = vb[itype][jtype];
        *gammaCur = gamma[itype][jtype];
        *thetaCur = theta[itype][jtype];
        const double Vstar = (*vbCur)/(std::pow(*R, 3));
        const double Sstar = (1+0.5*(*thetaCur))*(pow(Vstar,1/3.0) + 0.1*pow(Vstar,2.0/3.0)); // [Willett2000], equation (15), use the full-length e.g 2*Sc
        *critDist = Sstar*(*R);
        curCapModel = capModels[capillaryModelLoad[itype][jtype]-1];
          
        *firstTouchCap = 0;
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
      
      //fstat = DataFstat(Eigen::Vector3d P1, Eigen::Vector3d P2, int Id1, int Id2, Eigen::Vector3d Val, double VolWater, double DistCurr, double DistCrit)
    }

    void noCollision(ContactData& cdata, ForceData& i_forces, ForceData& j_forces){
      double * const history = &cdata.contact_history[history_offset];
      double * const firstTouch = &history[0];
      double * const touchFlag = &history[1];
      double * const firstTouchCap = &history[2];
      double * const critDist = &history[3];
      double * const R = &history[4];
      double * const vbCur = &history[5];
      double * const thetaCur = &history[6];
      double * const gammaCur = &history[7];
      
      if (*touchFlag and capillarFlag) {
        cdata.has_force_update = true;
        const Eigen::Vector3d dCur = Eigen::Vector3d(cdata.delta[0],cdata.delta[1],cdata.delta[2]);
        const Eigen::Vector3d dCurN = dCur.normalized();
        const double s = dCur.norm() - cdata.radsum;
        if(*critDist  > 0 and s < *critDist) {
          /*
          * Capillar model from Willet [Willett2000] (analytical solution), but
          * used also in the work of Herminghaus [Herminghaus2005]
          */
          const double Gamma = *gammaCur;
          const double Vb = *vbCur;
          const double Theta = *thetaCur;
          
          const double Fn = curCapModel(*R, Vb, Gamma, Theta, s); 
          
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
          *touchFlag = 0;
        }
      }
    }
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}
    DataFstat contactDataGet(CollisionData& cdata, ForceData&, ForceData&){
      //std::cout<<"EEEEEEEEEEEEEEEEE"<<std::endl;
      /*
      double * const history = &cdata.contact_history[history_offset];
      double * const firstTouch = &history[0];
      double * const touchFlag = &history[1];
      double * const firstTouchCap = &history[2];
      double * const critDist = &history[3];
      double * const R = &history[4];
      double * const vbCur = &history[5];
      double * const thetaCur = &history[6];
      double * const gammaCur = &history[7];
      */ 
    }

  protected:
    double ** k_n;
    double ** k_t;
    double ** gamma_n;
    double ** gamma_t;
    
    double ** gamma;
    double ** theta;
    double ** vb;
    double ** capillaryModelLoad;
    
    double ** e_n;
    double ** e_t;
    double ** t_c;
    double ** coeffFrict;
    
    bool explicitFlag;
    bool capillarFlag;
    
    int history_offset;
    
    std::vector<double (*)(const double &, const double &, const double &, const double &, const double &)> capModels;
    double (*curCapModel)(const double &, const double &, const double &, const double &, const double &);
    
  };
}
}


#endif // NORMAL_MODEL_HOOKE_VISCEL_H_
#endif
