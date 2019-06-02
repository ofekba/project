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

  allocate();
  if(update->ntimestep==0){
    printf("\n\n\n domain->box_exist=%d boxlo= %f %f %f boxhi= %f %f %f\n\n\n",domain->box_exist,domain->boxlo[0],domain->boxlo[1],domain->boxlo[2],domain->boxhi[0],domain->boxhi[1],domain->boxhi[2]);
  }
  box_x_len=abs(domain->boxlo[0])+abs(domain->boxhi[0]);
  box_y_len=abs(domain->boxlo[1])+abs(domain->boxhi[1]);
  box_z_len=abs(domain->boxlo[2])+abs(domain->boxhi[2]);
  printf("box bounderies %d %d %d\n",box_x_len,box_y_len,box_z_len);
  box_xhi=domain->boxhi[0];
  box_yhi=domain->boxhi[1];
  box_zhi=domain->boxhi[2];
  box_xlo=domain->boxlo[0];
  box_ylo=domain->boxlo[1];
  box_zlo=domain->boxlo[2];


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
  reaxc->set_extra_potential_parameters(); 
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

  FindNbr(lists);
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



void FixReaxCCheckFourset::FindNbr(struct _reax_list * /*lists*/)
{
 // printf("\n=========in FindNbr=========\n");//************
  int nlocal_tot = static_cast<int> (atom->nlocal);
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    //printf("\n==================nmax is=%d\n", nmax);//************
  }
  //printf("\n==================nlocal_tot=%d \n", nlocal_tot);
  

  for(int i=0; i<atom->nlocal; i++){
    int tag=atom->tag[i];
    tag_to_i[tag-1]=i;
  }

//NEIGH LIST BY CALCULATE DISTANCES

/*  double x0, x1, x2;
  double del0, del1, del2;
  double dist;
  for( int i = 0; i < atom->nlocal; ++i ){
    x0=atom->x[i][0];
    x1=atom->x[i][1];
    x2=atom->x[i][2];
    for( int j = 0; j < atom->nlocal; ++j ){
      del0=x0-atom->x[j][0];
      del1=x1-atom->x[j][1];
      del2=x2-atom->x[j][2];
      dist=del0*del0+del1*del1+del2*del2;
      dist=sqrt(dist);
      neigh_list[atom->tag[i]-1][atom->tag[j]-1]=dist;
    }
  }
  if(update->ntimestep<1){
    for(int i=8; i<10; i++){
      //nlocal_tot
        printf("___neigh of %d, num neigh:%d___\n", i+1,  numneigh[i]);
        for(int j=0; j<nlocal_tot; j++){
          printf("|id=%d, distance=%f",j+1, neigh_list[i][j]);
        }
        printf("|\n\n");
    }
  }*/


  /*for(int i=0; i<117; i++)
    for(int j=0; j<117; j++)
      neigh_list[i][j]=0.0;*/
  
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
      if(update->ntimestep>=11100 && update->ntimestep>=11105){
        if((atom_i->orig_id==28 && atom_i->orig_id==54) || (atom_i->orig_id==54 && atom_i->orig_id==28))
          printf("\n\n54-28 dist=%f\n", nbr_pj->d);
        if((atom_i->orig_id==54 && atom_i->orig_id==52) || (atom_i->orig_id==52 && atom_i->orig_id==54))
          printf("54-52 dist=%f\n", nbr_pj->d);
        if((atom_i->orig_id==52 && atom_i->orig_id==8) || (atom_i->orig_id==8 && atom_i->orig_id==52))
          printf("52-8 dist=%f\n", nbr_pj->d);
      }
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
  




    // PRINTS FOR DEBUGGING

    //print the neighbors list
  /* if(update->ntimestep<1){
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
            //printf("\nfirst cond OK\n");//**********
            for (c = 0; c < nlocal; c++) {
              c_tag = c+1; 
              i=tag_to_i[c_tag-1];
              if(i==-1)
                break;
              c_type = atom->type[i]; 
              if(c_type==THIRD_TYPE){
                if( 0.8<neigh_list[b][c] && neigh_list[b][c]<1.3 ){
                  //printf("\nsecond cond OK\n");//**********
                  for(d = 0; d < nlocal; d++) {
                    d_tag = d+1;
                    i=tag_to_i[d_tag-1];
                    if(i==-1)
                      break;
                    d_type = atom->type[i];
                    if(d_type==FORTH_TYPE){
                      if( 3.0<neigh_list[c][d] && neigh_list[c][d]<8.0 ){
                       // printf("\nthird cond OK\n");//**********
                        if( 0.9<neigh_list[d_tag-1][a_tag-1] && neigh_list[d_tag-1][a_tag-1]<2.2 ) {
                          //printf("\nfourth cond OK\n");//**********
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
                              //if(update->ntimestep>6000)
                                //printf("fourset- %d %d %d %d",a_tag,b_tag,c_tag,d_tag);
                            }
                            else if(neigh_list[a][b]+neigh_list[c][d]<n_9_min_oh_cn && c_tag==90){
                              n_9_min_oh_cn=neigh_list[a][b]+neigh_list[c][d];
                              fourset[1][0]=a_tag; //O
                              fourset[1][1]=b_tag; //H
                              fourset[1][2]=c_tag; //N
                              fourset[1][3]=d_tag; //C
                              num_fourset++;
                            }
                            //printf("**** success! ****\n") ; 

                            //random choice
                            /*fourset[num_fourset][0]=a_tag; //O
                            fourset[num_fourset][1]=b_tag; //H
                            fourset[num_fourset][2]=c_tag; //N
                            fourset[num_fourset][3]=d_tag; //C
                            num_fourset++;*/


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
 // printf("finish over the neigh list\n");

int apply_flag=0;
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

  //for level2 run
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
      printf("\nrand_num %d\n",rand_num);
      fprintf(fp,"# fourset O H N C at timestep " BIGINT_FORMAT " : ",update->ntimestep);
      fprintf(fp,"1/1- %d %d %d %d \n",fourset[0][0], fourset[0][1], fourset[0][2], fourset[0][3]);
    }

  }
    //FOR DEBUGGING
    /*if(apply_flag==1){
      double **x = atom->x;
      printf("\nOH=%f, CN=%f\n", neigh_list[fourset[0][0]-1][fourset[0][1]-1], neigh_list[fourset[0][2]-1][fourset[0][3]-1]);
      i=tag_to_i[fourset[0][0]-1];
      double xtmp = x[i][0];
      double ytmp = x[i][1];
      double ztmp = x[i][2];
      int j=tag_to_i[fourset[0][1]-1];
      double delx = xtmp - x[j][0];
      double dely = ytmp - x[j][1];
      double delz = ztmp - x[j][2];
      double rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
      double roh=sqrt(rsq);

      i=tag_to_i[fourset[0][2]-1];
       xtmp = x[i][0];
       ytmp = x[i][1];
       ztmp = x[i][2];
      j=tag_to_i[fourset[0][3]-1];
       delx = xtmp - x[j][0];
       dely = ytmp - x[j][1];
       delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz; //distance^2 between 2 atoms
      float rcn=sqrt(rsq);
      printf("\nOH=%f, CN=%f\n", roh, rcn);
    }*/
      
    //print the foursets
    //f(update->ntimestep<1){
    if(apply_flag==1){
      for(int i=0; i<num_fourset;i++){
        printf("fourset #%d is: ", i);
        for(int j=0; j<4;j++)
          printf("%d ",fourset[i][j]);
        printf("\n");
      }
      printf("\n");
    }
  

  //printf("\n\n==========finish checkFourset func=============\n");
}

/* ---------------------------------------------------------------------- */

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
  /*double x0, x1, x2;
  double j0, j1, j2;
  double del0, del1, del2;
  double dist;
  int num_of_epons=int(atom->nlocal/58.5);
  fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",update->ntimestep);
  for( int i = 0; i < atom->nlocal; ++i ){
    //level1-no if, level2-first, level2 noPBC-second line
    if(atom->tag[i]==8 ||atom->tag[i]==28 || atom->tag[i]==92 || atom->tag[i]==96 || atom->tag[i]==53 || atom->tag[i]==49 || atom->tag[i]==9|| atom->tag[i]==52|| atom->tag[i]==54|| atom->tag[i]==29|| atom->tag[i]==30|| atom->tag[i]==95|| atom->tag[i]==97){
    //if(atom->tag[i]==84 ||atom->tag[i]==90 || atom->tag[i]==60 || atom->tag[i]==76 || atom->tag[i]==59|| atom->tag[i]==88|| atom->tag[i]==102|| atom->tag[i]==103|| atom->tag[i]==99 || atom->tag[i]==91 || atom->tag[i]==94 ){

      fprintf(fp,"# atom %d type %d ",atom->tag[i], atom->type[i]);
      x0=atom->x[i][0];
      x1=atom->x[i][1];
      x2=atom->x[i][2];
      for( int j = 0; j < atom->nlocal; ++j ){
          j0=atom->x[j][0];
          j1=atom->x[j][1];
          j2=atom->x[j][2];

          del0=x0-j0;
          del1=x1-j1;
          del2=x2-j2;
          
         if( (abs(box_xhi-x0)<=3.2 && abs(box_xlo-j0)<=3.2) || (abs(box_xhi-j0)<=3.2 && abs(box_xlo-x0)<=3.2)){ 
            if(del0>0)
              del0-= box_x_len;
            else
              del0+= box_x_len;
          }
          if( (abs(box_yhi-x1)<=3.2 && abs(box_ylo-j1)<=3.2) || (abs(box_yhi-j1)<=3.2 && abs(box_ylo-x1)<=3.2)){
            if(del1>0)
              del1-= box_y_len;
            else
              del1+= box_y_len;
          }
          if( (abs(box_zhi-x2)<=3.2 && abs(box_zlo-j2)<=3.2) || (abs(box_zhi-j2)<=3.2 && abs(box_zlo-x2)<=3.2)){
            if(del2>0)
              del2-= box_z_len;
            else
              del2+= box_z_len;
          }

          dist=del0*del0+del1*del1+del2*del2;
          dist=sqrt(dist);
          fprintf(fp,"%d %f ",atom->tag[j], dist);
          //ofek
         // if(update->ntimestep>=10 && update->ntimestep<=20 && ((atom->tag[i]==53 && atom->tag[j]==70) || (atom->tag[i]==70 && atom->tag[j]==53)))
          //  printf("\n\t atom1 %d: x1 %f y1 %f z1 %f atom2 %d: x2 %f y2 %f z2 %f \ndelx %f dely %f delz %f dist %f\n",atom->tag[i], x0, x1, x2, atom->tag[j], atom->x[j][0], atom->x[j][1], atom->x[j][2], del0, del1, del2, dist);
          //if(j==i && update->ntimestep==50){
          //  printf(tag i=%d, )
         // }
      }
      fprintf(fp,"\n");
    }
  }
  fprintf(fp,"#\n");*/


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
   // if(fp!=NULL) fprintf(fp,"\n"); 
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


