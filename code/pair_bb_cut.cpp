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

pairBBcut::pairBBcut(LAMMPS *lmp) : Pair(lmp)
{
  allocate();
  int nmax = atom->nmax;
  f_fourset=NULL;
}

/* ---------------------------------------------------------------------- */

pairBBcut::~pairBBcut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(F1);
    memory->destroy(F2);
    memory->destroy(offset);
    memory->destroy(f_fourset);
  }
}

/* ---------------------------------------------------------------------- */

int pairBBcut::from_tag_to_i(tagint tag){
  if(tag<=0)
    return -1;
  //printf("----> in from_tag_to_i function <----\n") ; 
  //int len=sizeof(atom->tag)/sizeof(*atom->tag);
  for(int i=0; i<nmax; i++){
    if(int(atom->tag[i])==int(tag)){
        //printf("the tag is: %d the i is: %d\n", tag, i);
        //printf("----> out <----\n") ;
        return i;
    } 
  }
  
//printf("no i for that shitti tag");
  //printf("----> out no result with tag=%d<----\n", tag) ;
  return -1;

}





/* ---------------------------------------------------------------------- */

void pairBBcut::compute(int eflag, int vflag)
{
 /* for (int i = 0; i < nmax; i++) {
    for (int j = 0; j < 3; j++) {
      f_fourset[i][j] = 0;
    }
  }

  for(int k=0; k<num_fourset; k++){
    compute_pair(fourset[k][0], fourset[k][1]);
    compute_pair(fourset[k][0], fourset[k][3]);
    compute_pair(fourset[k][0], fourset[2][3]);
  }
  //reaxc.set_f_fourset(f_fourset);*/
 
}
/* ---------------------------------------------------------------------- */
//this method get:
//1. tag, type of atop i and tag, type of atom j
//2. the R(i,j)=the distance between them
//returns the calculated force.
double pairBBcut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq)
{
  double force;
  double r= sqrt(rsq)-wanted_dist[itype][jtype];
  double temp= -F2[itype][jtype] * r;
  force = -2 * F1[itype][jtype] * temp * exp(temp * r);
  return force;
}

/* ---------------------------------------------------------------------- */
void pairBBcut::compute_pair(int i_tag, int j_tag){
    
    int i,j,itype, jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz,fpair, rsq;
    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;

    i=from_tag_to_i(i_tag);

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    j=from_tag_to_i(j_tag);
    jtype = type[j];

    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
      
    rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
    fpair=single(i_tag, j_tag, itype, jtype, rsq);

    f_fourset[i_tag-1][0] += delx*fpair;
    f_fourset[i_tag-1][1] += dely*fpair;
    f_fourset[i_tag-1][2] += delz*fpair;
            
    f_fourset[j_tag-1][0] -= delx*fpair;
    f_fourset[j_tag-1][1] -= dely*fpair;
    f_fourset[j_tag-1][2] -= delz*fpair;

}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void pairBBcut::allocate()
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
  memory->create(f_fourset,nmax,3,"pair:f_fourset");//***************

  
  
}

/* ----------------------------------------------------------------------
                  -------->   TODO!!!!!!!    <--------  
   global settings
------------------------------------------------------------------------- */

void pairBBcut::settings(int narg, char **arg)
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

void pairBBcut::coeff(int narg, char **arg)
{
 
 //TODO:choose args for the coeff setting

 /* if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();*/
  
  for(int i=0; i<atom->ntypes+1; i++){
    for(int j=0; j<atom->ntypes+1; j++){
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

void pairBBcut::init_style()
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

void *pairBBcut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"F1") == 0 || strcmp(str,"F1") == 0) return (void *) F1;
  if (strcmp(str,"F2") == 0 || strcmp(str,"F2") == 0) return (void *) F2;
  return NULL;
}
