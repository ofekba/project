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
#include <cstdlib>
#include "fix.h"
#include "pointers.h"

namespace LAMMPS_NS {

class FixReaxCCheckFourset : public Fix {
 public:
  FixReaxCCheckFourset(class LAMMPS *, int, char **); //constructor
  virtual ~FixReaxCCheckFourset(); //destructor
  int setmask();
  virtual void init();
  void setup(int);
  void end_of_step();
  int MAX_NUM_FOURSETS;

 protected:
  int me, nprocs, nmax, ntypes, maxsize;
  int **fourset; //list of fourset to apply the potential on
  int **o_c_pair_tags; //the legal O-C atom pairs pattern to apply the extra potential on
  int *n_tags; //the pattern of N atom's tags
  int num_fourset; //0 if the list is empty. else, number of foursets in the list
  FILE *fp; //for dists file that follow the distances between atoms.
  char *fp_suffix; //dists file suffix
  int set_mol_pattern();//function to set the pattern of o_c_pair_tags and n_tags
  int nevery_dists_follow; //nevery for the dists documentation
  int nevery_file_dists; //seperate the dists file into many files with this timesteps nevery
  int nevery_cond_check; // every this many steps look for legal foursets

  void allocate(); //alocate memory
  void destroy(); //free memory
  virtual void Output_ReaxC_Bonds(bigint);
  void FindNbr(struct _reax_list*);
  int nint(const double &);
  virtual double memory_usage();

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxC *reaxc;
  class NeighList *list;
};
}

#endif
#endif