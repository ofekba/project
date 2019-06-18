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
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::FixReaxCCheckFourset(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //checking for correct usage with Fix ReaxC Check_Fourset command in the "in" file
  if (narg != 5) error->all(FLERR,"Illegal fix reax/c/checkFourset command");
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;
  nmax = atom->nmax;
  
  nevery = force->inumeric(FLERR,arg[3]);

  if (nevery <= 0 )
    error->all(FLERR,"Illegal fix reax/c/checkFourset command, illigal nevery");//**********

//for dists fp that follow the distance between atoms
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
  fourset = NULL;

  allocate();
}

/* ---------------------------------------------------------------------- */

FixReaxCCheckFourset::~FixReaxCCheckFourset()
{
  MPI_Comm_rank(world,&me);
  destroy();
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCCheckFourset::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::init()
{
  //reaxc singelton instance
  reaxc = (PairReaxC *) force->pair_match("reax/c",0);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/checkFourset without "
                                "pair_style reax/c, reax/c/kk, or reax/c/omp");
  reaxc->set_extra_potential_parameters(); //set the parameter of the extra potential from the user input file
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::end_of_step()
{
  Output_ReaxC_Bonds(update->ntimestep);
  followDistFunc();
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::Output_ReaxC_Bonds(bigint /*ntimestep*/)
{
  int i, j;
  //reset the neigh struct 
  for (i = 0; i < atom->nlocal; i++) {
    for (j = 0; j < atom->nlocal; j++) {
      neigh_list[i][j]=-1;
    }
  }
  //check for distances
  FindNbr(lists);
  //edit distance file
  checkForFoursets();

}

/* ---------------------------------------------------------------------- */



void FixReaxCCheckFourset::FindNbr(struct _reax_list * /*lists*/)
{
  int nlocal_tot = static_cast<int> (atom->nlocal); //number of atoms in the system
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
  }
  
//update tag to i array (to convert tag to i index in the atom list)
  for(int i=0; i<atom->nlocal; i++){
    int tag=atom->tag[i];
    tag_to_i[tag-1]=i;
  }
  
  //NEIGH LIST BY FAR NEIGH LIST STRUCT
  int i, j, pj;
  int start_i, end_i;
  int type_i, type_j;
  int btop_i;
  reax_list *far_nbrs;
  far_neighbor_data *nbr_pj;
  reax_atom *atom_i, *atom_j;

  far_nbrs = (reaxc->lists) + FAR_NBRS;

  for( i = 0; i < reaxc->system->N; ++i ) {
    atom_i = &(reaxc->system->my_atoms[i]);
    type_i  = atom_i->type;
    if (type_i < 0) continue;
    
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);

    double **x = atom->x;
    double xtmp,ytmp,ztmp,delx,dely,delz,rsq,rij;
    int jjj, iii;
    
    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      j = nbr_pj->nbr;
      atom_j = &(reaxc->system->my_atoms[j]);
  
      if (nbr_pj->d <= 8.0){
        //add to neigh list
        iii=tag_to_i[atom_i->orig_id-1];
        xtmp = x[iii][0];
        ytmp = x[iii][1];
        ztmp = x[iii][2];
        jjj=tag_to_i[atom_j->orig_id-1];
        delx = xtmp - x[jjj][0];
        dely = ytmp - x[jjj][1];
        delz = ztmp - x[jjj][2];
        rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
        rij=sqrt(rsq);
        if(rij<8.0){
          neigh_list[atom_i->orig_id-1][atom_j->orig_id-1]=nbr_pj->d;
          neigh_list[atom_j->orig_id-1][atom_i->orig_id-1]=nbr_pj->d;
      }

        
      }
    }
  }
  

}

/* ---------------------------------------------------------------------- */
/* search for foursets that meets the paper conditions for legal fourset
to apply the extra potential on, and the condition on O-C pairs
(C atom that bonded only to one O atom) */
void FixReaxCCheckFourset::checkForFoursets(){
  
  //printf("\n============in checkFourset===============\n");
  int nlocal = atom->nlocal;

  //reset the fourset array
  num_fourset=0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 4; j++) {
      fourset[i][j] = 0;
    }
  }

  //mine
  tagint a_tag, b_tag, c_tag, d_tag;
  int a_type, b_type, c_type, d_type;
  int a, ai, b, c, d;
  int i; //the index of an atom in the atom->tag array

  //const
  int const FIRST_TYPE=3; //O TYPE NUM
  int const SECOND_TYPE=2; // H TYPE  NUM
  int const THIRD_TYPE=4;// N TYPE  NUM
  int const FORTH_TYPE=1;// C TYPE  NUM
  
  //choose fourset to apply the potential on with the smallest OH CN distances
  double n_8_min_oh_cn=17;
  double n_9_min_oh_cn=17;
  double min_oh_cn=17;


  for (ai = 0; ai < nlocal; ai++) {
    a_tag = atom->tag[ai];
    a_type = atom->type[ai];
    a=a_tag-1;
    if(a_type==FIRST_TYPE){
      for (b = 0; b < nlocal; b++) {
        b_tag = b+1;
        i=tag_to_i[b_tag-1];
        if(i==-1)
          break;
        b_type = atom->type[i]; 
        if(b_type==SECOND_TYPE){
          if( 1.3<neigh_list[a][b] && neigh_list[a-1][b]<8.0 ){
            //found legal O-H
            for (c = 0; c < nlocal; c++) {
              c_tag = c+1; 
              i=tag_to_i[c_tag-1];
              if(i==-1)
                break;
              c_type = atom->type[i]; 
              if(c_type==THIRD_TYPE){
                if( 0.8<neigh_list[b][c] && neigh_list[b][c]<1.3 ){
                  //found legal N-H
                  for(d = 0; d < nlocal; d++) {
                    d_tag = d+1;
                    i=tag_to_i[d_tag-1];
                    if(i==-1)
                      break;
                    d_type = atom->type[i];
                    if(d_type==FORTH_TYPE){
                      if( 3.0<neigh_list[c][d] && neigh_list[c][d]<8.0 ){
                       //found legal N-C
                        if( 0.9<neigh_list[d_tag-1][a_tag-1] && neigh_list[d_tag-1][a_tag-1]<2.2 ) {
                          //found legal C-O
                          
                          //code for level1 run
                          /*if(neigh_list[a][b]+neigh_list[c][d]<n_8_min_oh_cn){
                            n_8_min_oh_cn=neigh_list[a][b]+neigh_list[c][d];
                            num_fourset++;

                            fourset[num_fourset][0]=fourset[0][0]; //O
                            fourset[num_fourset][1]=fourset[0][1]; //H
                            fourset[num_fourset][2]=fourset[0][2]; //N
                            fourset[num_fourset][3]=fourset[0][3]; //C

                            fourset[0][0]=a_tag; //O
                            fourset[0][1]=b_tag; //H
                            fourset[0][2]=c_tag; //N
                            fourset[0][3]=d_tag; //C
                          }*/

                          //code for level2 run N=8/9
                          //if(a_tag==d_tag+4 || a_tag==d_tag+2){
                          
                          //for noPBC run N=84,90.
                          if( (a_tag==60 && d_tag==59) || (a_tag==116 && d_tag==113) || (a_tag==99 && d_tag==94) || (a_tag==103 && d_tag==102) ){
                            if(neigh_list[a][b]+neigh_list[c][d]<n_8_min_oh_cn && c_tag==84){
                              n_8_min_oh_cn=neigh_list[a][b]+neigh_list[c][d];
                              fourset[0][0]=a_tag; //O
                              fourset[0][1]=b_tag; //H
                              fourset[0][2]=c_tag; //N
                              fourset[0][3]=d_tag; //C
                              num_fourset++;
                            }
                            else if(neigh_list[a][b]+neigh_list[c][d]<n_9_min_oh_cn && c_tag==90){
                              n_9_min_oh_cn=neigh_list[a][b]+neigh_list[c][d];
                              fourset[1][0]=a_tag; //O
                              fourset[1][1]=b_tag; //H
                              fourset[1][2]=c_tag; //N
                              fourset[1][3]=d_tag; //C
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
          }
        }
      }
    }
  }


int apply_flag=0; //turn to 1 if the reaxc apply the extra potential on the chosen fourset. else, 0
  //operate the extra potential
  if(num_fourset>0){
    
    //for level1 run
    /*apply_flag=reaxc->set_fourset(fourset, 1);
    if(apply_flag==1){
      printf("list of foursets:\n");
      for(int i=0; i<num_fourset;i++){
          printf("fourset #%d is: ", i);
          for(int j=0; j<4;j++)
            printf("%d ",fourset[i][j]);
          printf("\n");
      }
      printf("\n");
      fprintf(fp,"# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
      fprintf(fp,"1/1- %d %d %d %d \n",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
    }
  }*/

    //extra for noPBC run
  if(update->ntimestep<1002 || update->ntimestep>11000){
      if(n_8_min_oh_cn<17){
        apply_flag=reaxc->set_fourset(fourset, 1);
        if(apply_flag==1){
          fprintf(fp,"# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
          fprintf(fp,"1/1- %d %d %d %d \n",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
        }
        return;
      }
      else
        return;
    }

  //for any level2 run
    if(n_8_min_oh_cn<17 && n_9_min_oh_cn<17){
      if(n_8_min_oh_cn>n_9_min_oh_cn){
        fourset[0][0]=fourset[1][0]; //O
        fourset[0][1]=fourset[1][1]; //H
        fourset[0][2]=fourset[1][2]; //N
        fourset[0][3]=fourset[1][3]; //C
      }
    }
    else{
      if(n_9_min_oh_cn<17){
        fourset[0][0]=fourset[1][0]; //O
        fourset[0][1]=fourset[1][1]; //H
        fourset[0][2]=fourset[1][2]; //N
        fourset[0][3]=fourset[1][3]; //C
      }
    }
    apply_flag=reaxc->set_fourset(fourset, 1);
      if(apply_flag==1){
        fprintf(fp,"# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
        fprintf(fp,"1/1- %d %d %d %d \n",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
      }
    
    //random choice
    int rand_num;
    /*if(num_fourset>1){
       rand_num= int(rand() % num_fourset + 1) - 1;
       if(update->ntimestep>11000 && update->ntimestep<11100) rand_num=1;
      int temp;
      //if(update->ntimestep>12000 && update->ntimestep<13000) rand_num=0;
      temp=fourset[0][0];
      fourset[0][0] = fourset[rand_num][0];
      fourset[rand_num][0]=temp;
      temp=fourset[0][1];
      fourset[0][1] = fourset[rand_num][1];
      fourset[rand_num][1]=temp;
      temp=fourset[0][2];
      fourset[0][2] = fourset[rand_num][2];
      fourset[rand_num][2]=temp;
      temp=fourset[0][3];
      fourset[0][3] = fourset[rand_num][3];
      fourset[rand_num][3]=temp;
    }*/
    /*if(fourset[0][0]==103 && fourset[0][2]==84 && fourset[0][1]==88 && fourset[0][3]==102){
      if(update->ntimestep%100 == 0)
      printf("\nnum foursets %d\n",num_fourset);
      return;
    }*/
    apply_flag=reaxc->set_fourset(fourset, 1);
    if(apply_flag==1){
      fprintf(fp,"# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
      fprintf(fp,"1/1- %d %d %d %d \n",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
    }

  }
      
    //OPTIONAL: print the foursets he found.
    if(apply_flag==1){
      for(int i=0; i<num_fourset;i++){
        printf("fourset #%d is: ", i);
        for(int j=0; j<4;j++)
          printf("%d ",fourset[i][j]);
        printf("\n");
      }
      printf("\n");
    }
  
}

/* ---------------------------------------------------------------------- */
/* document distances between each two atoms to the dists file */
void FixReaxCCheckFourset::followDistFunc()
{
  int _nevery=10;
  if(update->ntimestep==0){
    fprintf(fp,"# totalTimesteps " BIGINT_FORMAT " \n",update->laststep);
    fprintf(fp,"# totalAtomNum %d \n",atom->nlocal);
    fprintf(fp,"# fix_nevery %d \n",_nevery);
  }
  if(update->ntimestep%_nevery!=0)
    return;

//NEIGH LIST BY FAR NEIGH LIST STRUCT
  if(fp!=NULL) fprintf(fp,"# Timestep " BIGINT_FORMAT " ",update->ntimestep);
  int pi1;
  int i1, i2;
  int start_i1, end_i1;
  int type_i1, type_i2, tag_i1, tag_i2;
  reax_list *far_nbrs, *bond_nbrs;
  far_neighbor_data *nbr_p_far;
  bond_data *nbr_p_bond;
  reax_atom *atom_i1, *atom_i2;

  far_nbrs = (reaxc->lists) + FAR_NBRS;
  bond_nbrs = (reaxc->lists) + BONDS;

  for( i1 = 0; i1 < reaxc->system->N; ++i1 ) {
    atom_i1 = &(reaxc->system->my_atoms[i1]);
    type_i1  = atom_i1->type;
    tag_i1 = atom_i1->orig_id;
      
    //level2
    if(tag_i1==84 ||tag_i1==90 || tag_i1==60 || tag_i1==76 || tag_i1==59|| tag_i1==88|| tag_i1==102|| tag_i1==103|| tag_i1==99 || tag_i1==91 || tag_i1==94 ){
    //if(tag_i1==8 ||tag_i1==28 || tag_i1==92 || tag_i1==96 || tag_i1==53 || tag_i1==49 || tag_i1==9|| tag_i1==52|| tag_i1==54|| tag_i1==29|| tag_i1==30|| tag_i1==95|| tag_i1==97){
      if(fp!=NULL){
          fprintf(fp,"\n# atom %d type %d ",tag_i1, type_i1+1);
      }

      start_i1 = Start_Index(i1, far_nbrs);
      end_i1   = End_Index(i1, far_nbrs);
      for( pi1 = start_i1; pi1 < end_i1; ++pi1 ) {
        nbr_p_far = &( far_nbrs->select.far_nbr_list[pi1] );
        i2 = nbr_p_far->nbr;
        atom_i2 = &(reaxc->system->my_atoms[i2]);
        type_i2= atom_i2->type;
        tag_i2 = atom_i2->orig_id;
        if(fp!=NULL){
            fprintf(fp,"%d %f ",tag_i2, nbr_p_far->d);
        }
      }
      
      start_i1 = Start_Index(i1, bond_nbrs);
      end_i1   = End_Index(i1, bond_nbrs);
      for( pi1 = start_i1; pi1 < end_i1; ++pi1 ) {
        nbr_p_bond = &( bond_nbrs->select.bond_list[pi1] );
        i2 = nbr_p_bond->nbr;
        atom_i2 = &(reaxc->system->my_atoms[i2]);
        type_i2= atom_i2->type;
        tag_i2 = atom_i2->orig_id;
        if(fp!=NULL){
            fprintf(fp,"%d %f ",tag_i2, nbr_p_bond->d);
        }
      }
    }
  }
  if(fp!=NULL) fprintf(fp,"\n#\n");
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
  memory->destroy(fourset);
  memory->destroy(tag_to_i);
}

/* ---------------------------------------------------------------------- */

void FixReaxCCheckFourset::allocate()
{
  memory->create(fourset,20,4,"reax/c/checkFourset:fourset");//***************
  memory->create(neigh_list,atom->nlocal,atom->nlocal,"reax/c/checkFourset:neigh_list");
  memory->create(tag_to_i,atom->nlocal,"reax/c/checkFourset:tag_to_i");

}

/* ---------------------------------------------------------------------- */

double FixReaxCCheckFourset::memory_usage()
{
  double bytes;
  bytes += atom->nlocal*sizeof(int);//tag_to_i
  bytes += 1.0*atom->nlocal*4*sizeof(int);//fourset
  bytes += 1.0*atom->nlocal*atom->nlocal*sizeof(double);//neigh_list
  return bytes;
}


