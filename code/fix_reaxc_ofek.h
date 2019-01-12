/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(reax/c/ofek,FixReaxCOfek)

#else

#ifndef LMP_FIX_REAXC_OFEK_H
#define LMP_FIX_REAXC_OFEK_H

#include <cstdio>
#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixReaxCOfek : public Fix {
 public:
  FixReaxCOfek(class LAMMPS *, int, char **);
  virtual ~FixReaxCOfek();
  int setmask();
  virtual void init();
  void setup(int);
  void end_of_step();

 protected:
  int me, nprocs, nmax, ntypes, maxsize;
  int *numneigh; // list of numneigh for each atom
  tagint **neighid; //list of neigh per atom by ID
  double **neigh_d; //distance between neighbors
  int **fourset; //list of fourset to appky tha potential on
  int num_fourset; //0 if the list is empty. else, number of fourset
  int local_tot; //num of atoms

  void allocate();
  void destroy();
  virtual void Output_ReaxC_Bonds(bigint);
  void FindNbr(struct _reax_list*, int &);
  int nint(const double &);
  virtual double memory_usage();
  void OfekFunc(); //*****my func*****
  int from_tag_to_i(tagint tag);//*****my func*****

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxC *reaxc;
  class NeighList *list;
};
}

#endif
#endif
