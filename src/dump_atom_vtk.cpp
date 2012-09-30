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

#include "string.h"
#include "dump_atom_vtk.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpATOMVTK::DumpATOMVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump filename");

  size_one = 5;
  sort_flag = 1;
  sortcol = 0;

  char *str = (char *) "%d %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::write_header(bigint n)
{
}

/* ---------------------------------------------------------------------- */

int DumpATOMVTK::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::pack(int *ids)
{
  int n = 0;
  
  vtkSmartPointer<vtkPoints>  spheresPos = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();
  
  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("radii");
  
  vtkSmartPointer<vtkDoubleArray> spheresMass = vtkSmartPointer<vtkDoubleArray>::New();
  spheresMass->SetNumberOfComponents(1);
  spheresMass->SetName("mass");
  
  vtkSmartPointer<vtkIntArray> spheresId = vtkSmartPointer<vtkIntArray>::New();
  spheresId->SetNumberOfComponents(1);
  spheresId->SetName("id");
  
  vtkSmartPointer<vtkIntArray> spheresType = vtkSmartPointer<vtkIntArray>::New();
  spheresType->SetNumberOfComponents(1);
  spheresType->SetName("type");
  
  vtkSmartPointer<vtkDoubleArray> spheresVelL = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelL->SetNumberOfComponents(3);
  spheresVelL->SetName("velocity_lin");
  
  vtkSmartPointer<vtkDoubleArray> spheresVelA = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelA->SetNumberOfComponents(3);
  spheresVelA->SetName("velocity_ang");
  
  vtkSmartPointer<vtkDoubleArray> spheresForce = vtkSmartPointer<vtkDoubleArray>::New();
  spheresForce->SetNumberOfComponents(3);
  spheresForce->SetName("force");
  
  
  

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **o = atom->omega;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (ids) ids[n++] = tag[i];
      
      vtkIdType pid[1];
      pid[0] = spheresPos->InsertNextPoint(x[i][0], x[i][1], x[i][2]);
      radii->InsertNextValue(atom->radius[i]);
      
      double vv[3] = { (double) v[i][0], (double) v[i][1], (double) v[i][2] };
      spheresVelL->InsertNextTupleValue(vv);
      
      double oo[3] = { (double) o[i][0], (double) o[i][1], (double) o[i][2] };
      spheresVelA->InsertNextTupleValue(oo);

      double ff[3] = { (double) f[i][0], (double) f[i][1], (double) f[i][2] };
      spheresForce->InsertNextTupleValue(ff);
      
      if (atom->mass) {
        spheresMass->InsertNextValue(atom->mass[i]);
      } else {
        spheresMass->InsertNextValue(0.0);
      }
      
      spheresId->InsertNextValue(tag[i]);
      spheresType->InsertNextValue(type[i]);
      
      spheresCells->InsertNextCell(1,pid);
    }
  }
  
  spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  spheresUg->SetPoints(spheresPos);
  spheresUg->SetCells(VTK_VERTEX, spheresCells);
  spheresUg->GetPointData()->AddArray(radii);
  spheresUg->GetPointData()->AddArray(spheresId);
  spheresUg->GetPointData()->AddArray(spheresType);
  spheresUg->GetPointData()->AddArray(spheresMass);
  spheresUg->GetPointData()->AddArray(spheresVelL);
  spheresUg->GetPointData()->AddArray(spheresVelA);
  spheresUg->GetPointData()->AddArray(spheresForce);
  
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::write_data(int n, double *mybuf)
{
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  
  writer->SetInput(spheresUg);
  
  
  //HACK!!! The stream should be set instead of temporary file!!!
  //===================================
  writer->SetFileName("tmp.vtu");
  writer->Write();
  
  std::string line;
  std::ifstream myfile ("tmp.vtu");
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      getline (myfile,line);
      char *a=new char[line.size()+1];
      a[line.size()]=0;
      memcpy(a,line.c_str(),line.size());
      fprintf(fp,"%s\n",a);
    }
    myfile.close();
  } else std::cout << "Unable to open file"; 
  remove( "tmp.vtu" );
  //===================================
  
}
