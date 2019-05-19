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
  //printf("\n*********FixReaxCCheckFourset:\tin constructor***********\n");
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

  fourset = NULL;

  allocate();
  //printf("\n*********FixReaxCCheckFourset:\tout constructor***********\n");

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
  //printf("\n****out setmask*****\n");
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
  //printf("\n*********FixReaxCCheckFourset:\tin init***********\n");
  reaxc = (PairReaxC *) force->pair_match("reax/c",0);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/checkFourset without "
                                "pair_style reax/c, reax/c/kk, or reax/c/omp");
  //printf("\n*********FixReaxCCheckFourset:\tout allocate***********\n");
}

/* ---------------------------------------------------------------------- */
//function that create a file that following the distance between all of the atoms.
void FixReaxCCheckFourset::end_of_step()
{
  //printf("\n****in end_of_step*****\n");
  Output_ReaxC_Bonds(update->ntimestep);
  //followDistFunc();
  if (me == 0) fflush(fp);
  //printf("\n****out end_of_step*****\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::Output_ReaxC_Bonds(bigint /*ntimestep*/)
{
  FindNbr(lists);
}

/* ---------------------------------------------------------------------- */


/*
   - list->inum  = the length of the neighborlists list.
   - list->ilist  = list list of "i" atoms for which neighbor lists exist.
   - list->numneigh = the length of each these neigbor list.
>  - list->firstneigh = the pointer to the list of neighbors of "i".
*/



void FixReaxCCheckFourset::FindNbr(struct _reax_list * /*lists*/)
{
  //printf("\n=========in FindNbr=========\n");//************
  int nlocal_tot = static_cast<int> (atom->nlocal);
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
  }
  int _nevery=10;
  int num_of_epons=int(atom->nlocal/58.5);
  if(fp!=NULL){
    if(update->ntimestep==0){
      fprintf(fp,"# totalTimesteps " BIGINT_FORMAT " \n",update->laststep);
      fprintf(fp,"# totalAtomNum %d \n",atom->nlocal);
      fprintf(fp,"# fix_nevery %d",_nevery);
    }
    if(update->ntimestep%_nevery==0)
      fprintf(fp,"\n# Timestep " BIGINT_FORMAT ,update->ntimestep);
  }

  if(update->ntimestep%_nevery!=0)
    return;

  const int TYPE_C = 0;
  const int TYPE_H = 1;
  const int TYPE_O = 2;
  const int TYPE_N = 3;
  
  //NEIGH LIST BY FAR NEIGH LIST STRUCT
  int pi1, pi2, pi3, pi4;
  int _o, _h, _n, _c;
  int start_o, end_o, start_h, end_h, start_n, end_n, start_c, end_c, start_12, end_12;
  int type_i1, type_i2, type_i3, type_i4, tag_i1, tag_i2;
  reax_list *far_nbrs, *bond_nbrs;
  far_neighbor_data *nbr_p_oh, *nbr_p_nc;
  bond_data *nbr_p_hn, *nbr_p_co, *nbr_p_12;
  reax_atom *atom_i1, *atom_i2, *atom_i3, *atom_i4, *atom_i5;

  far_nbrs = (reaxc->lists) + FAR_NBRS;
  bond_nbrs = (reaxc->lists) + BONDS;
  num_fourset = 0;

  for(int nn=0; nn<atom->nlocal; nn++){
    fourset[nn][0]=0;
    fourset[nn][1]=0;
    fourset[nn][2]=0;
    fourset[nn][3]=0;
  }

  for( _o = 0; _o < reaxc->system->N; ++_o ) {
    atom_i1 = &(reaxc->system->my_atoms[_o]);
    type_i1  = atom_i1->type;
    tag_i1 = atom_i1->orig_id;
    if (type_i1 != TYPE_O) continue;
    start_o = Start_Index(_o, far_nbrs);
    end_o   = End_Index(_o, far_nbrs);
    start_12 = Start_Index(_o, bond_nbrs);
    end_12   = End_Index(_o, bond_nbrs);

    if(fp!=NULL){
      if((tag_i1-8)%43==0 ||(tag_i1-16)%43==0  || (tag_i1-22)%43==0 ==60 || (tag_i1-18)%43==0  || (tag_i1-23)%43==0 || (tag_i1-21)%43==0 || (tag_i1-15)%43==0 || (tag_i1-19)%43==0 || (tag_i1-(num_of_epons*43)-8)%31==0 || (tag_i1-(num_of_epons*43)-9)%31==0 ){
        fprintf(fp,"\n# atom %d type %d ",tag_i1, type_i1);
      }
    }

    for( pi1 = start_12; pi1 < end_12; ++pi1 ) {
      nbr_p_12 = &( bond_nbrs->select.bond_list[pi1] );
      _h = nbr_p_12->nbr;
      atom_i2 = &(reaxc->system->my_atoms[_h]);
      type_i2= atom_i2->type;
      tag_i2 = atom_i2->orig_id;
      if(fp!=NULL && tag_i2>0){
        if((tag_i1-8)%43==0 ||(tag_i1-16)%43==0  || (tag_i1-22)%43==0 ==60 || (tag_i1-18)%43==0  || (tag_i1-23)%43==0 || (tag_i1-21)%43==0 || (tag_i1-15)%43==0 || (tag_i1-19)%43==0 || (tag_i1-(num_of_epons*43)-8)%31==0 || (tag_i1-(num_of_epons*43)-9)%31==0 ){
          fprintf(fp,"%d %f ",tag_i2, nbr_p_12->d);
        }
      }
    }

    for( pi1 = start_o; pi1 < end_o; ++pi1 ) {
      nbr_p_oh = &( far_nbrs->select.far_nbr_list[pi1] );
      _h = nbr_p_oh->nbr;
      atom_i2 = &(reaxc->system->my_atoms[_h]);
      type_i2= atom_i2->type;
      tag_i2 = atom_i2->orig_id;

      if(fp!=NULL){
        if((tag_i1-8)%43==0 ||(tag_i1-16)%43==0  || (tag_i1-22)%43==0 ==60 || (tag_i1-18)%43==0  || (tag_i1-23)%43==0 || (tag_i1-21)%43==0 || (tag_i1-15)%43==0 || (tag_i1-19)%43==0 || (tag_i1-(num_of_epons*43)-8)%31==0 || (tag_i1-(num_of_epons*43)-9)%31==0 ){
          fprintf(fp,"%d %f ",tag_i2, nbr_p_oh->d);
        }
      }

      //printf("o %d type %d atom %d type %d\n", atom_i1->orig_id, type_i1, atom_i2->orig_id, type_i2);
      if(type_i2 != TYPE_H) continue;
      if (1.3 <= nbr_p_oh->d && nbr_p_oh->d <= 8.0 ){
        //printf("\n\n*****cond 1 OK*****\n\n");
        start_h = Start_Index(_h, bond_nbrs);
        end_h = End_Index(_h, bond_nbrs);
        for( pi2 = start_h; pi2 < end_h; ++pi2 ){
          nbr_p_hn = &( bond_nbrs->select.bond_list[pi2] );
          _n=nbr_p_hn->nbr;
          atom_i3 = &(reaxc->system->my_atoms[_n]);
          type_i3= atom_i3->type;
         //printf("h %d type %d atom %d type %d\n", atom_i2->orig_id,type_i2, atom_i3->orig_id, type_i3);
          if(type_i3 != TYPE_N) continue;
          if (0.8 <= nbr_p_hn->d && nbr_p_hn->d <= 1.3 ){
            //printf("\n\n*****cond 2 OK doki*****\n\n");
            start_n = Start_Index(_n, far_nbrs);
            end_n = End_Index(_n, far_nbrs);
            for( pi3 = start_n; pi3 < end_n; ++pi3 ){
              nbr_p_nc = &( far_nbrs->select.far_nbr_list[pi3] );
              _c=nbr_p_nc->nbr;
              atom_i4 = &(reaxc->system->my_atoms[_c]);
              type_i4= atom_i4->type;
              if(type_i4 != TYPE_C) continue;
              if(3.0 <= nbr_p_nc->d && nbr_p_nc->d <= 8.0 ){
                //printf("*****cond 3 OK doki*****\n");

                //options for o-c
                //o= 22+43x c=18+43x || o= 8+43x c=16+43x || o= 23+43x c=21+43x || o= 15+43x c=19+43x
                int _re=atom_i1->orig_id%43;
                int _x=int( (atom_i1->orig_id - _re) / 43 );
                int _optional_c_tag=0;
                switch(_re) {
                  case 22: _optional_c_tag=18+43*_x;
                    break;
                  case 8: _optional_c_tag=16+43*_x;
                    break;
                  case 23: _optional_c_tag=21+43*_x;
                    break;
                  case 15: _optional_c_tag=19+43*_x;
                    break;
                }
                if(atom_i4->orig_id != _optional_c_tag) break;

                start_c = Start_Index(_c, bond_nbrs);
                end_c = End_Index(_c, bond_nbrs);
                
                for( pi4 = start_c; pi4 < end_c; ++pi4 ){
                  nbr_p_co = &( bond_nbrs->select.bond_list[pi4] );
                  atom_i5=&(reaxc->system->my_atoms[nbr_p_co->nbr]);
                  if(nbr_p_co->nbr != _o) continue;
                  if (0.9 <= nbr_p_co->d && nbr_p_co->d <= 2.2 ){
                    //printf("\n\n*****cond 4 OK *****\n\n");
                    //printf("\n\n*~~~~*\t\tsuccess !!!!!!!\t\t *~~~~*\n\n");
                    fourset[num_fourset][0] = atom_i1->orig_id;
                    fourset[num_fourset][1] = atom_i2->orig_id;
                    fourset[num_fourset][2] = atom_i3->orig_id;
                    fourset[num_fourset][3] = atom_i4->orig_id;
                    num_fourset++;

                  }
                }
              }
            }
          }
        }
      }
    }
  }
//OFEK
  /*if(update->ntimestep%100==0){
    printf("\n\nat timestep %d\n",update->ntimestep);
    for(int nn=0; nn<num_fourset; nn++)
      printf("fourset #%d: %d %d %d %d\n",nn, fourset[nn][0], fourset[nn][1], fourset[nn][2], fourset[nn][3]);
    printf("\n");
  }*/

    if(num_fourset!=0){
      int rand_num= int(rand() % num_fourset + 1) - 1;
      //printf("~~~~~~rand forset=%d~~~~~~\n", rand_num);
      fourset[0][0] = fourset[rand_num][0];
      fourset[0][1] = fourset[rand_num][1];
      fourset[0][2] = fourset[rand_num][2];
      fourset[0][3] = fourset[rand_num][3];

      //OFEK
      //TODO: to operate the potentil on 2 fourset togther
      
      /*int rand_flag=0;
      int sec_rand_num;
      int is_already_found[num_fourset];
      for(int nn=0; nn<num_fourset; nn++)
        is_already_found[nn]=0;
      is_already_found[rand_num-1]=1;
      while(rand_flag<num_fourset){
        sec_rand_num= int(rand() % num_fourset + 1) - 1;
        rand_flag++;
        if(is_already_found[sec_rand_num-1]==0){
          if()

        }
        
      }*/

      int apply_flag = reaxc->set_fourset(fourset, 1);
      //OFEK
      if(apply_flag==1)
        printf("\nstart operate the potential\n");
    }
  
  

      //printf("\n=============finish FindNbr=========\n");
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
  memory->destroy(fourset);
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::allocate()
{
  //printf("\n*********FixReaxCCheckFourset:\tin allocate***********\n");
  memory->create(fourset,atom->nlocal,4,"reax/c/checkFourset:fourset");//***************
  //printf("\n*********FixReaxCCheckFourset:\tout allocate***********\n");
}

/* ---------------------------------------------------------------------- */

double FixReaxCCheckFourset::memory_usage()
{
  double bytes;
  //bytes = 3.0*nmax*sizeof(double);//??
  bytes += 1.0*atom->nlocal*4*sizeof(int);//fourset

  return bytes;
}


