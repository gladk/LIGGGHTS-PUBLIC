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

#ifdef DUMP_CLASS

DumpStyle(atom/vtk,DumpATOMVTK)

#else

#ifndef LMP_DUMP_ATOM_VTK_H
#define LMP_DUMP_ATOM_VTK_H

#include "dump.h"

#include<vtkCellArray.h>
#include<vtkFloatArray.h>
#include<vtkDoubleArray.h>
#include<vtkIntArray.h>
#include<vtkPoints.h>
#include<vtkPointData.h>
#include<vtkCellData.h>
#include<vtkSmartPointer.h>
#include<vtkUnstructuredGrid.h>
#include<vtkPolyData.h>
#include<vtkXMLUnstructuredGridWriter.h>
#include<vtkXMLPolyDataWriter.h>
#include<vtkZLibDataCompressor.h>
#include<vtkTriangle.h>
#include<vtkLine.h>
#include<vtkQuad.h>

namespace LAMMPS_NS {

class DumpATOMVTK : public Dump {
 public:
  DumpATOMVTK(class LAMMPS *, int, char**);
  ~DumpATOMVTK() {}

 private:
  void init_style();
  void write_header(bigint);
  int count();
  void pack(int *);
  void write_data(int, double *);
  
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg;
};

}

#endif
#endif
