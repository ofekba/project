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

FixStyle(reax/c/checkFourset,FixReaxCCheckFourset)

#else

#ifndef LMP_FIX_REAXC_CHECK_FOURSET_H
#define LMP_FIX_REAXC_CHECK_FOURSET_H

#include <cstdio>
#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixReaxCCheckFourset : public Fix {
 public:
  FixReaxCCheckFourset(class LAMMPS *, int, char **);
  virtual ~FixReaxCCheckFourset();
  int setmask();
  virtual void init();
  void setup(int);
  void end_of_step();

 protected:
  int me, nprocs, nmax, ntypes, maxsize;
  int **fourset; //list of fourset to appky tha potential on
  int num_fourset; //0 if the list is empty. else, number of fourset
  double **neigh_list;
  int *tag_to_i;

  void allocate();
  void destroy();
  virtual void Output_ReaxC_Bonds(bigint);
  void FindNbr(struct _reax_list*);
  int nint(const double &);
  virtual double memory_usage();
  void checkForFoursets(); //*****my func*****
  int from_tag_to_i(tagint tag);//*****my func*****
  FILE *fp;
  void followDistFunc();

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxC *reaxc;
  class NeighList *list;
};
}

#endif
#endif
