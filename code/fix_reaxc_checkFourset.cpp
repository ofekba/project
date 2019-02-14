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
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include "fix_ave_atom.h"
#include "fix_reaxc_checkFourset.h"
#include "atom.h"
#include "update.h"
#include "pair_reaxc.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::FixReaxCCheckFourset(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix reax/c/checkFourset command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;
  nmax = atom->nmax;
  //printf("\n==================nmax=%d\n", nmax);//************
  
  nevery = force->inumeric(FLERR,arg[3]);

  if (nevery <= 0 )
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, illigal nevery");//**********

//for fp
  if (me == 0) {
      char *suffix = strrchr(arg[4],'.');
      if (suffix && strcmp(suffix,".gz") == 0) {
  #ifdef LAMMPS_GZIP
        char gzip[128];
        snprintf(gzip,128,"gzip -6 > %s",arg[4]);
  #ifdef _WIN32
        fp = _popen(gzip,"wb");
  #else
        fp = popen(gzip,"w");
  #endif
  #else
        error->one(FLERR,"Cannot open gzipped file");
  #endif
      } else fp = fopen(arg[4],"w");

      if (fp == NULL) {
        char str[128];
        snprintf(str,128,"Cannot open fix reax/c/bonds file %s",arg[4]);
        error->one(FLERR,str);
      }
    }

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for fix reax/c bonds");

  neigh_list=NULL;
  tag_to_i=NULL;
  numneigh = NULL;
  fourset = NULL;
  num_fourset = 0;

  allocate();

  for(int i=0; i<atom->nlocal; i++){
    int tag=atom->tag[i];
    tag_to_i[tag]=i;
  }
}

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::~FixReaxCCheckFourset()
{
  //printf("\n****in destructor*****\n");
  MPI_Comm_rank(world,&me);
  destroy();
  if (me == 0) fclose(fp);
  //printf("\n*-*-*-*finish the fix of ofki ofkilish :)*-*-*-*\n");
}

/* ---------------------------------------------------------------------- */

int FixReaxCCheckFourset::setmask()
{
  //printf("\n*****in setmask*****\n");
  int mask = 0;
  mask |= END_OF_STEP;
 // printf("\n****out setmask*****\n");
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::setup(int /*vflag*/)
{
  //printf("\n****in setup*****\n");
  end_of_step();
  //printf("\n****out setup*****\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::init()
{
  reaxc = (PairReaxC *) force->pair_match("reax/c",0);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/checkFourset without "
                                "pair_style reax/c, reax/c/kk, or reax/c/omp");

}

/* ---------------------------------------------------------------------- */
//function that create a file that following the distance between all of the atoms.
void FixReaxCCheckFourset::end_of_step()
{
  //printf("\n****in end_of_step*****\n");
  Output_ReaxC_Bonds(update->ntimestep);
  followDistFunc();
  if (me == 0) fflush(fp);
  //printf("\n****out end_of_step*****\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::Output_ReaxC_Bonds(bigint /*ntimestep*/)
{
  //printf("\n****in Output_ReaxC_Bonds*****\n");
  int i, j;
  int numbonds;

  //printf("\n****** 1 ******\n");
  /*if (atom->nmax > nmax) {
    destroy();
    nmax = atom->nmax;
     printf("\n==================nmax num 2=%d\n", nmax);//************
    allocate();
  }*/
  for (i = 0; i < atom->nlocal; i++) {
    numneigh[i] = 0;
    for (j = 0; j < atom->nlocal; j++) {
      neigh_list[i][j]=-1;
    }
  }
 // printf("\n****** 2 ******\n");
  numbonds = 0;
  num_fourset=0;
  for (i = 0; i < 20; i++) {
    for (j = 0; j < 4; j++) {
      fourset[i][j] = 0;
    }
  }
//printf("\n****** 3 ******\n");
  FindNbr(lists, numbonds);
  //printf("\n****** 4 ******\n");
  checkForFoursets();

  //printf("\n=======finish Output_ReaxC_Bonds=====\n");//************Nbr

}

/* ---------------------------------------------------------------------- */


/*
   - list->inum  = the length of the neighborlists list.
   - list->ilist  = list list of "i" atoms for which neighbor lists exist.
   - list->numneigh = the length of each these neigbor list.
>  - list->firstneigh = the pointer to the list of neighbors of "i".
*/



void FixReaxCCheckFourset::FindNbr(struct _reax_list * /*lists*/, int &numbonds)
{
 // printf("\n=========in FindNbr=========\n");//************
  int nlocal_tot = static_cast<int> (atom->nlocal);
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    //printf("\n==================nmax is=%d\n", nmax);//************
  }
  //printf("\n==================nlocal_tot=%d \n", nlocal_tot);
  
  int i, j, pj;
  int tag_i, tag_j;
  far_neighbor_data *nbr_pj;
  int start_i, end_i;
  int type_i, type_j;
  reax_atom *atom_i, *atom_j;
  reax_list *far_nbrs = (reaxc->lists)+ FAR_NBRS;
  double cutoff = reaxc->control->bond_cut;


  for (i = 0; i < reaxc->system->N; i++) {
    //printf("\n=============i=%d=========\n",i);
    atom_i = &(reaxc->system->my_atoms[i]);
    type_i  = atom_i->type;
    tag_i=(int)atom_i->orig_id;
    if(tag_i <= 0||tag_i > nlocal_tot)
      continue;
    start_i = Start_Index(i, far_nbrs);
    end_i  = End_Index(i, far_nbrs);
    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      if (nbr_pj->d <= cutoff) {
        j = nbr_pj->nbr;
        atom_j = &(reaxc->system->my_atoms[j]);
        if(atom_j->orig_id > 0 && atom_j->orig_id<=nlocal_tot){
          tag_j=atom_j->orig_id;
          if(neigh_list[tag_i-1][tag_j-1] == -1){
            numneigh[tag_i-1]++;
            numneigh[tag_j-1]++;
          }
          neigh_list[tag_i-1][tag_j-1]=nbr_pj->d;
          neigh_list[tag_j-1][tag_i-1]=nbr_pj->d;
         }
      }
    }
  }

//print the neighbors list
  /*if(update->ntimestep<99){
    for(i=0; i<nlocal_tot; i++){
        printf("___neigh of %d, num neigh:%d___\n", i+1,  numneigh[i]);
        for(j=0; j<nlocal_tot; j++){
          printf("|id=%d, distance=%f",j+1, neigh_list[i][j]);
        }
        printf("|\n\n");
    }
  }*/
    //  printf("\n=============finish FindNbr=========\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::checkForFoursets(){

  int nlocal = atom->nlocal;
  //printf("\n============in checkFourset===============\n");
  
  //mine
  tagint a_tag, b_tag, c_tag, d_tag;
  int a_type, b_type, c_type, d_type;
  int a, b,c,d, a_numNbr, b_numNbr, c_numNbr, d_numNbr;
  int i; //the index of an atom in the atom->tag array

  //const
  int const FIRST_TYPE=3; //O TYPE NUM
  int const SECOND_TYPE=2; // H TYPE  NUM
  int const THIRD_TYPE=4;// N TYPE  NUM
  int const FORTH_TYPE=1;// C TYPE  NUM
  for (a = 0; a < nlocal; a++) {
    a_tag = atom->tag[a];
    a_type = atom->type[a];
    int success=0;
    if(a_type==FIRST_TYPE){
      a_numNbr = nint(numneigh[a_tag-1]);
      for (b = 0; b < a_numNbr; b++) {
       /* if(success==1)
          break;*/
        b_tag = b+1;
        i=tag_to_i[b_tag];
        if(i==-1)
          break;
        b_type = atom->type[i]; 
        if(b_type==SECOND_TYPE){
          if( 1.3<neigh_list[a][b] && neigh_list[a][b]<8.0 ){
            b_numNbr=nint(numneigh[b_tag-1]); 
            //printf("\nfirst cond OK\n");//**********

            for (c = 0; c < b_numNbr; c++) {
              /*if(success==1)
                break;*/
              c_tag = c+1; 
              i=tag_to_i[c_tag];
              if(i==-1)
                break;
              c_type = atom->type[i]; 
              if(c_type==THIRD_TYPE){
                if( 0.8<neigh_list[b][c] && neigh_list[b][c]<1.3 ){
                  printf("\nsecond cond OK\n");//**********

                  c_numNbr=nint(numneigh[c_tag-1]);
                  for(d = 0; d < c_numNbr; d++) {
                    /*if(success==1)
                      break;*/
                    d_tag = d+1;
                    i=tag_to_i[d_tag];
                    if(i==-1)
                      break;
                    d_type = atom->type[i];
                    
                    if(d_type==FORTH_TYPE){
                      if( 3.0<neigh_list[c][d] && neigh_list[c][d]<8.0 ){
                        printf("\nthird cond OK\n");//**********
                        
                        if( 0.9<neigh_list[d_tag-1][a_tag-1] && neigh_list[d_tag-1][a_tag-1]<2.2 ) {
                          //printf("\nfourth cond OK\n");//**********
                          
                          num_fourset++;
                          fourset[num_fourset-1][0]=a_tag; //O
                          fourset[num_fourset-1][1]=b_tag; //H
                          fourset[num_fourset-1][2]=c_tag; //N
                          fourset[num_fourset-1][3]=d_tag; //C
                          //printf("**** success! ****\n") ; 
                          success=1;
                          //break;
                        }
                      }
                    } 
                  }
                }
              }
            }
          }
        }
      }
    }
  }
 // printf("finish over the neigh list\n");
  if(num_fourset!=0){
    //reaxc->set_fourset(fourset, num_fourset);
    for(int i=0; i<num_fourset;i++){
      printf("fourset #%d is: ", num_fourset);
      for(int j=0; j<4;j++)
        printf("%d ",fourset[i][j]);
      printf("\n");
    }
  }
  //printf("\n\n==========finish checkFourset func=============\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::followDistFunc()
{
  double x0, x1, x2;
  double del0, del1, del2;
  double dist;
  fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",update->ntimestep);
  for( int i = 0; i < atom->nlocal; ++i ){
    fprintf(fp,"# atom %d type %d ",atom->tag[i], atom->type[i]);
    x0=atom->x[i][0];
    x1=atom->x[i][1];
    x2=atom->x[i][2];
    for( int j = 0; j < atom->nlocal; ++j ){
      
        del0=x0-atom->x[j][0];
        del1=x1-atom->x[j][1];
        del2=x2-atom->x[j][2];
        dist=del0*del0+del1*del1+del2*del2;
        dist=sqrt(dist);
        fprintf(fp,"%d %f ",atom->tag[j], dist);
        /*if(j==i && update->ntimestep==50){
          printf(tag i=%d, )
        }*/
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"#\n");
}

/* ---------------------------------------------------------------------- */


int FixReaxCCheckFourset::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::destroy()
{
  memory->destroy(neigh_list);
  memory->destroy(numneigh);
  memory->destroy(fourset);
  memory->destroy(tag_to_i);
  

}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::allocate()
{
  memory->create(fourset,20,4,"reax/c/checkFourset:fourset");//***************
  memory->create(neigh_list,atom->nlocal,atom->nlocal,"reax/c/checkFourset:neigh_list");
  memory->create(numneigh,atom->nlocal,"reax/c/checkFourset:numneigh");
  memory->create(tag_to_i,atom->nlocal,"reax/c/checkFourset:tag_to_i");
  
}

/* ---------------------------------------------------------------------- */

double FixReaxCCheckFourset::memory_usage()
{
  double bytes;
  //bytes = 3.0*nmax*sizeof(double);//??
  bytes = atom->nlocal*sizeof(int);//numneigh
  bytes += atom->nlocal*sizeof(int);//tag_to_i
  bytes += 1.0*atom->nlocal*4*sizeof(int);//fourset
  bytes += 1.0*atom->nlocal*atom->nlocal*sizeof(double);//neigh_list

  return bytes;
}


