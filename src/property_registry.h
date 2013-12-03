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
#ifndef PROPERTY_REGISTRY_H_
#define PROPERTY_REGISTRY_H_

#include "pointers.h"
#include "lammps.h"
#include "fix_property_global.h"
#include <map>
#include <set>
#include <string>
#include "mech_param_gran.h"
#include "error.h"
#include "modify.h"

using namespace std;
using namespace LAMMPS_NS;

class PropertyRegistry;

template<typename T>
class Property {
public:
  Property() : data(0) {}
  virtual ~Property(){}

  void connect(T & variable) {
    listeners.insert(&variable);
    variable = data;
  }

  void disconnect(T & variable) {
    typename set<T*>::iterator it = listeners.find(&variable);
    listeners.erase(it);
    variable = NULL;
  }

  void updateAll() {
    for(typename set<T*>::iterator it = listeners.begin(); it != listeners.end(); ++it) {
      *(*it) = data;
    }
  }

  T data;
  set<T*> listeners;
};

typedef Property<double> ScalarProperty;

class VectorProperty : public Property<double*>
{
public:
  VectorProperty(const int N){
    data = new double[N];
    for(int col = 0; col < N; col++) {
      data[col] = 0.0;
    }
  }

  virtual ~VectorProperty(){
    delete [] data;
  }
};

class MatrixProperty : public Property<double**>
{
public:
  MatrixProperty(const int N, const int M){
    double * array = new double[N*M];
    data = new double*[N];

    for(int row = 0; row < N; row++) {
      data[row] = &array[row*M];

      for(int col = 0; col < M; col++) {
        data[row][col] = 0.0;
      }
    }
  }

  virtual ~MatrixProperty() {
    delete [] data[0];
    delete [] data;
  }
};

// -------------------------------------------------------------------

typedef ScalarProperty* (*ScalarPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);
typedef VectorProperty* (*VectorPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);
typedef MatrixProperty* (*MatrixPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);

class PropertyRegistry : protected Pointers {
  MechParamGran mpg;

public:
  PropertyRegistry(LAMMPS* lmp) : Pointers(lmp), mpg(lmp)
  {
  }

  ~PropertyRegistry()
  {
  }

  int max_type() {
    return mpg.max_type();
  }

  LAMMPS * getLAMMPS() {
    return lmp;
  }

  FixPropertyGlobal* getGlobalProperty(const char *varname,const char *style,const char *svmstyle,int len1,int len2,const char *caller){
    return static_cast<FixPropertyGlobal*>(modify->find_fix_property(varname, style, svmstyle, len1, len2, caller));
  }

  ScalarProperty * getScalarProperty(string varname) {
    if(scalars.find(varname) == scalars.end()) {
      if(scalar_creators.find(varname) != scalar_creators.end()) {
        scalars[varname] = (*scalar_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        error->message(FLERR, "unknown scalar property");
      }
    }
    return scalars[varname];
  }

  VectorProperty * getVectorProperty(string varname) {
    if(vectors.find(varname) == vectors.end()) {
      if(vector_creators.find(varname) != vector_creators.end()) {
        vectors[varname] = (*vector_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        error->message(FLERR, "unknown vector property");
      }
    }
    return vectors[varname];
  }

  MatrixProperty * getMatrixProperty(string varname) {
    if(matrices.find(varname) == matrices.end()) {
      if(matrix_creators.find(varname) != matrix_creators.end()) {
        matrices[varname] = (*matrix_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        error->message(FLERR, "unknown matrix property");
      }
    }
    return matrices[varname];
  }

  void registerProperty(string varname, ScalarPropertyCreator creator, bool sanity_checks = false) {
    if(scalar_creators.find(varname) == scalar_creators.end()) {
      scalar_creators[varname] = creator;
      use_sanity_checks[varname] = sanity_checks;
    } else if(scalar_creators[varname] != creator) {
      error->message(FLERR, "property with the same name, but different implementation registered");
    }
  }

  void registerProperty(string varname, VectorPropertyCreator creator, bool sanity_checks = false) {
    if(vector_creators.find(varname) == vector_creators.end()) {
      vector_creators[varname] = creator;
      use_sanity_checks[varname] = sanity_checks;
    } else if(vector_creators[varname] != creator) {
      error->message(FLERR, "property with the same name, but different implementation registered");
    }
  }

  void registerProperty(string varname, MatrixPropertyCreator creator, bool sanity_checks = false) {
    if(matrix_creators.find(varname) == matrix_creators.end()) {
      matrix_creators[varname] = creator;
      use_sanity_checks[varname] = sanity_checks;
    } else if(matrix_creators[varname] != creator) {
      error->message(FLERR, "property with the same name, but different implementation registered");
    }
  }

  void connect(string varname, double ** & variable) {
    if(matrices.find(varname) == matrices.end()) {
      if(matrix_creators.find(varname) != matrix_creators.end()) {
        matrices[varname] = (*matrix_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        // ERROR unknown property
      }
    }
    matrices[varname]->connect(variable);
  }

  void connect(string varname, double * & variable) {
    if(vectors.find(varname) == vectors.end()) {
      if(vector_creators.find(varname) != vector_creators.end()) {
        vectors[varname] = (*vector_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        // ERROR unknown property
      }
    }
    vectors[varname]->connect(variable);
  }

  void connect(string varname, double & variable) {
    if(scalars.find(varname) == scalars.end()) {
      if(scalar_creators.find(varname) != scalar_creators.end()) {
        scalars[varname] = (*scalar_creators[varname])(*this, "", use_sanity_checks[varname]);
      } else {
        // ERROR unknown property
      }
    }
    scalars[varname]->connect(variable);
  }

  void init() {
    scalars.clear();
    vectors.clear();
    matrices.clear();
  }

private:
  map<string, ScalarPropertyCreator> scalar_creators;
  map<string, VectorPropertyCreator> vector_creators;
  map<string, MatrixPropertyCreator> matrix_creators;

  map<string, ScalarProperty*> scalars;
  map<string, VectorProperty*> vectors;
  map<string, MatrixProperty*> matrices;

  map<string, bool> use_sanity_checks;
};

#endif /* PROPERTY_REGISTRY_H_ */
