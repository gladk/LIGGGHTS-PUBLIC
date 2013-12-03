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
#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_EPSD,epsd,2)
#else
#ifndef ROLLING_MODEL_EPSD_H_
#define ROLLING_MODEL_EPSD_H_
#include "contact_models.h"
#include <algorithm>
#include "math.h"
#include "domain.h"
#include "math_extra_liggghts.h"

namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<typename Style>
  class RollingModel<ROLLING_EPSD, Style> : protected Pointers
  {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    RollingModel(class LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp), coeffRollFrict(NULL), coeffRollVisc(NULL)
    {
      history_offset = hsetup->add_value("r_torquex_old", "1");
      hsetup->add_value("r_torquey_old", "1");
      hsetup->add_value("r_torquez_old", "1");
      STATIC_ASSERT(Style::TANGENTIAL == TANGENTIAL_HISTORY);
      printf("EPSD loaded\n");
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.registerProperty("coeffRollVisc", &MODEL_PARAMS::createCoeffRollVisc);
      registry.connect("coeffRollFrict", coeffRollFrict);
      registry.connect("coeffRollVisc", coeffRollVisc);
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) 
    {
      double r_torque[3], wr_roll[3];
      vectorZeroize3D(r_torque);

      if(cdata.touch) *cdata.touch |= TOUCH_ROLLING_MODEL;
      double * const c_history = &cdata.contact_history[history_offset]; // requires Style::TANGENTIAL == TANGENTIAL_HISTORY

      const double rmu= coeffRollFrict[cdata.itype][cdata.jtype];

      if(cdata.is_wall) {
        double dr_torque[3],wr_n[3],wr_t[3];

        const int itype = cdata.itype;
        const double enx = cdata.en[0];
        const double eny = cdata.en[1];
        const double enz = cdata.en[2];
        const double wr1 = cdata.wr1;
        const double wr2 = cdata.wr2;
        const double wr3 = cdata.wr3;
        const double radius = cdata.radi;
        double dt = update->dt; 

        // remove normal (torsion) part of relative rotation
        // use only tangential parts for rolling torque
        double wr_dot_delta = wr1*enx+ wr2*eny + wr3*enz;
        wr_n[0] = enx * wr_dot_delta;
        wr_n[1] = eny * wr_dot_delta;
        wr_n[2] = enz * wr_dot_delta;
        wr_t[0] = wr1 - wr_n[0];
        wr_t[1] = wr2 - wr_n[1];
        wr_t[2] = wr3 - wr_n[2];

        // spring
        const double kr = 2.25*cdata.kn*rmu*rmu*radius*radius; 
        // TODO ROLLING != 2 => kr = kt*radius*radius;

        dr_torque[0] = kr * wr_t[0] * dt;
        dr_torque[1] = kr * wr_t[1] * dt;
        dr_torque[2] = kr * wr_t[2] * dt;

        r_torque[0] = c_history[0] + dr_torque[0];
        r_torque[1] = c_history[1] + dr_torque[1];
        r_torque[2] = c_history[2] + dr_torque[2];

        // limit max. torque
        const double r_torque_mag = vectorMag3D(r_torque);
        const double r_torque_max = fabs(cdata.Fn)*radius*rmu;
        if(r_torque_mag > r_torque_max)
        {
          //printf("[%d] %e > %e\n", update->ntimestep, r_torque_mag, r_torque_max);
          const double factor = r_torque_max / r_torque_mag;

          r_torque[0] *= factor;
          r_torque[1] *= factor;
          r_torque[2] *= factor;

          // save rolling torque due to spring
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];

          // no damping / no dashpot in case of full mobilisation rolling angle
          // r_coef = 0.0;
        } else {
          //printf("[%d] %e <= %e\n", update->ntimestep, r_torque_mag, r_torque_max);
          // save rolling torque due to spring before adding damping torque
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];

          // dashpot
          double r_inertia;
          if (domain->dimension == 2) r_inertia = 1.5*cdata.mi*radius*radius;
          else  r_inertia = 1.4*cdata.mi*radius*radius;

          const double r_coef = coeffRollVisc[itype][cdata.jtype] * 2 * sqrt(r_inertia*kr);

          // add damping torque
          r_torque[0] += r_coef*wr_t[0];
          r_torque[1] += r_coef*wr_t[1];
          r_torque[2] += r_coef*wr_t[2];
        }
      } else {
        double wr_roll_n[3],wr_roll_t[3];

        const int i = cdata.i;
        const int j = cdata.j;
        const int itype = cdata.itype;
        const int jtype = cdata.jtype;
        const double enx = cdata.en[0];
        const double eny = cdata.en[1];
        const double enz = cdata.en[2];

        const double radi = cdata.radi;
        const double radj = cdata.radj;
        const double * const * const omega = atom->omega;

        // relative rotational velocity
        vectorSubtract3D(omega[i],omega[j],wr_roll);

        // remove normal (torsion) part of relative rotation
        // use only tangential parts for rolling torque
        double wr_dot_delta = wr_roll[0]*enx+ wr_roll[1]*eny + wr_roll[2]*enz;
        wr_roll_n[0] = enx * wr_dot_delta;
        wr_roll_n[1] = eny * wr_dot_delta;
        wr_roll_n[2] = enz * wr_dot_delta;
        vectorSubtract3D(wr_roll,wr_roll_n,wr_roll_t);

        // spring
        const double reff=radi*radj/(radi+radj);
        const double kr = 2.25*cdata.kn*rmu*rmu*reff*reff; 

        double dr_torque[3];
        vectorScalarMult3D(wr_roll_t,update->dt*kr,dr_torque); 

        r_torque[0] = c_history[0] + dr_torque[0];
        r_torque[1] = c_history[1] + dr_torque[1];
        r_torque[2] = c_history[2] + dr_torque[2];

        // limit max. torque
        const double r_torque_mag = vectorMag3D(r_torque);
        const double r_torque_max = fabs(cdata.Fn)*reff*rmu;
        if(r_torque_mag > r_torque_max)
        {
          const double factor = r_torque_max / r_torque_mag;

          r_torque[0] *= factor;
          r_torque[1] *= factor;
          r_torque[2] *= factor;

          // save rolling torque due to spring
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];
        } else {
          // save rolling torque due to spring before adding damping torque
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];

          // dashpot
          const double r_inertia_red_i = cdata.mi*radi*radi;
          const double r_inertia_red_j = cdata.mj*radj*radj;
          double r_inertia;
          if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
          else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);

          const double r_coef = coeffRollVisc[itype][jtype] * 2 * sqrt(r_inertia*kr);

          // add damping torque
          r_torque[0] += r_coef*wr_roll_t[0];
          r_torque[1] += r_coef*wr_roll_t[1];
          r_torque[2] += r_coef*wr_roll_t[2];
        }
      }

      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];
      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];
    }

    void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      if(cdata.touch) *cdata.touch &= ~TOUCH_ROLLING_MODEL;
      double * const c_history = &cdata.contact_history[history_offset];
      c_history[0] = 0.0; // this is the r_torque_old
      c_history[1] = 0.0; // this is the r_torque_old
      c_history[2] = 0.0; // this is the r_torque_old
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;
    double ** coeffRollVisc;
    int history_offset;
  };
}
#endif // ROLLING_MODEL_EPSD_H_
#endif
