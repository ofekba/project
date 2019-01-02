/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_bb_cut.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

pair_bb_cut::pair_bb_cut(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

pair_bb_cut::~pair_bb_cut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(F1);
    memory->destroy(F2);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void pair_bb_cut::compute(double **Foursets, int num_foursets)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq, temp, r12, rij;
  int *ilist,*jlist,*numneigh,**firstneigh;


  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for(int k=0; ik<num_foursets; k++){
    for(int q=0; q<2; q++)
    i = Foursets[k][q];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    j=i = Foursets[k][q+1];
    jtype = type[j];

    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
    
    rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
    rij=sqrt(rsq);
    r12=wanted_dist[itype][jtype];

    double r= rij-r12;
    double temp= -r*F2[itype][jtype];
    fpair = 2 * F1[itype][jtype] * temp * exp(temp * r)

    f[i][0] += delx*fpair;
    f[i][1] += dely*fpair;
    f[i][2] += delz*fpair;
          
    f[j][0] -= delx*fpair;
    f[j][1] -= dely*fpair;
    f[j][2] -= delz*fpair;

  }
  
  
  // loop over neighbors of my atoms
  /*for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
      rij=sqrt(rsq);
      r12=wanted_dist[itype][jtype];
 
      double r= sqr_dist-r12;
      double temp= -r*F2[itype][jtype];
      fpair = 2 * F1[itype][jtype] * temp * exp(temp * r)

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
          
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;
    
    }
  }*/
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void pair_bb_cut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  //needed???
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  //mine
  memory->create(wanted_dist,n+1,n+1,"pair:wanted_dist");
  memory->create(F1,n+1,n+1,"pair:F1");
  memory->create(F2,n+1,n+1,"pair:F2");

  
  
}

/* ----------------------------------------------------------------------
                  -------->   TODO!!!!!!!    <--------  
   global settings
------------------------------------------------------------------------- */

void pair_bb_cut::settings(int narg, char **arg)
{
  /*if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }*/
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void pair_bb_cut::coeff(int narg, char **arg)
{
 
 //TODO:choose args for the coeff setting

 /* if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();*/
  
  for(int i=0; i<atom->ntypes; i++){
    for(int j=0; j<atom->ntypes; j++){
      F1[i][j]=F2[i][j]=0;
      wanted_dist[i][j]=0;
    }
  }
  /*  C_type=1 H_type=2 O_type=3 N_type=4   */

  //O-C (3-1)
  F1[3][1]=F1[1][3]=50;
  F2[3][1]=F2[1][3]=0.5;
  wanted_dist[3][1]=wanted_dist[1][3]=1.95;

  //O-H (3-2)
  F1[3][2]=F1[2][3]=250;
  F2[3][2]=F2[2][3]=0.75;
  wanted_dist[3][2]=wanted_dist[2][3]=1.1;

  //N-C (4-1)
  F1[4][1]=F1[1][4]=300;
  F2[4][1]=F2[1][4]=0.75;
  wanted_dist[4][1]=wanted_dist[1][4]=1.5;
  
}

/* ----------------------------------------------------------------------
                -------->   TODO!!!!!!!    <--------  
   init specific to this pair style
------------------------------------------------------------------------- */

void pair_bb_cut::init_style()
{
  // request regular or rRESPA neighbor list

  int irequest;
  int respa = 0;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  irequest = neighbor->request(this,instance_me);

  if (respa >= 1) {
    neighbor->requests[irequest]->respaouter = 1;
    neighbor->requests[irequest]->respainner = 1;
  }
  if (respa == 2) neighbor->requests[irequest]->respamiddle = 1;

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */
double pair_bb_cut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq)
{
  double force,philj;

  double r= sqrt(rsq)-wanted_dist[itype][jtype]
  double temp= -F2[itype][jtype] * r;
  force = -2 * F1[itype][jtype] * temp * exp(temp * r)
  return force;
}

/* ---------------------------------------------------------------------- */

void *pair_bb_cut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"F1") == 0 || strcmp(str,"F1") == 0) return (void *) F1;
  if (strcmp(str,"F2") == 0 || strcmp(str,"F2") == 0) return (void *) F2;
  return NULL;
}
