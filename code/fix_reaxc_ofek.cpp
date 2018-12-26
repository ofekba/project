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
#include "fix_reaxc_ofek.h"
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

FixReaxCOfek::FixReaxCOfek(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix reax/c/ofek command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;
  nmax = atom->nmax;
  printf("\n==================nmax=%d\n", nmax);//************
  
  nevery = force->inumeric(FLERR,arg[3]);

  if (nevery <= 0 )
    error->all(FLERR,"Illegal fix reax/c/ofek command, illigal nevery");//**********

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
      snprintf(str,128,"Cannot open fix reax/c/ofek file %s",arg[4]);
      error->one(FLERR,str);
    }
  }

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for fix reax/c bonds");

  neigh_d= NULL;
  neighid = NULL;
  numneigh = NULL;
  local_tot = static_cast<int> (atom->natoms);

  allocate();
}

/* ---------------------------------------------------------------------- */

FixReaxCOfek::~FixReaxCOfek()
{
  MPI_Comm_rank(world,&me);

  destroy();

  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCOfek::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::init()
{
  reaxc = (PairReaxC *) force->pair_match("reax/c",0);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/ofek without "
                                "pair_style reax/c, reax/c/kk, or reax/c/omp");

}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::end_of_step()
{
  Output_ReaxC_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::Output_ReaxC_Bonds(bigint /*ntimestep*/, FILE * /*fp*/)

{
  int i, j;
  int nbuf, nbuf_local;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);

  /*if (atom->nmax > nmax) {
    destroy();
    nmax = atom->nmax;
     printf("\n==================nmax num 2=%d\n", nmax);//************
    allocate();
  }*/
  for (i = 0; i < nmax; i++) {
    numneigh[i] = 0;
    for (j = 0; j < nmax; j++) {
      neighid[i][j] = 0;
      neigh_d[i][j]=0.0;//**************
    }
  }
  numbonds = 0;

  FindNbr(lists, numbonds);

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);

  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memory->create(buf,nbuf,"reax/c/ofek:buf");
  for (i = 0; i < nbuf; i ++) buf[i] = 0.0;

  // Pass information to buffer
  PassBuffer(buf, nbuf_local);

  // Receive information from buffer for output
  RecvBuffer(buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memory->destroy(buf);

}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::FindNbr(struct _reax_list * /*lists*/, int &numbonds)
{
     
     int nlocal_tot = static_cast<int> (atom->natoms);

     printf("\n==================nlocal_tot=%d\n", nlocal_tot);
  int i, j, pj, nj;
  far_neighbor_data *nbr_pj;
  int start_i, end_i;
  int type_i, type_j;
  reax_atom *atom_i, *atom_j;
  reax_list *far_nbrs = (reaxc->lists)+ FAR_NBRS;;
  double cutoff = reaxc->control->bond_cut;

 /* int numNeighPerAtom[local_tot+1];
  for(i=0; i<local_tot+1; i++)
    numNeighPerAtom[i]=0;*/

printf("\n==================\n");
for(i = 0; i < 100; i++){
    atom_i = &(reaxc->system->my_atoms[i]);
    int tag_i=(int)atom_i->orig_id;
    type_i  = atom_i->type;
    printf("\ni=%d, orig_id=%d, imprt_id=%d, type=%d",i,tag_i, atom_i->imprt_id, type_i);
}
printf("\n==================\n");

  for (i = 0; i < nlocal_tot; i++) {
    printf("\n=============i=%d=========\n",i);
    atom_i = &(reaxc->system->my_atoms[i]);
    int tag_i=(int)atom_i->orig_id;
    type_i  = atom_i->type;
    start_i = Start_Index(i, far_nbrs);
    end_i  = End_Index(i, far_nbrs);
    nj=0;
    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
      if (nbr_pj->d <= cutoff) {
        j = nbr_pj->nbr;
        atom_j = &(reaxc->system->my_atoms[j]);
        type_j = atom_j->type;
        neighid[tag_i-1][nj] = atom_j->orig_id;
         //printf("\nd is=%f in [%d] [%d]", nbr_pj->d, tag_i ,nj);//***************
         neigh_d[tag_i-1][nj]=nbr_pj->d;

        nj++;
      }
    }
    numneigh[tag_i-1]=nj;
    //numNeighPerAtom[tag_i]=nj;
  }
      printf("\n=============finish=========\n");
    for(i=0; i<nlocal_tot; i++){
        printf("neigh of %d:", i+1);
        for(j=0; j<numneigh[i]; j++){
            printf("|id=%d, distance=%f",int(neighid[i][j]), neigh_d[i][j]);
        }
        printf("|\n\n");
    }
}


/* ---------------------------------------------------------------------- */

void FixReaxCOfek::PassBuffer(double *buf, int &nbuf_local)
{
  int i, j, k, numNbrs;
  int nlocal = atom->nlocal;

  j = 2;
  buf[0] = nlocal;
  for (i = 0; i < nlocal; i++) {
    buf[j-1] = atom->tag[i];
    buf[j+0] = atom->type[i];
    buf[j+2] = reaxc->workspace->nlp[i];
    buf[j+3] = atom->q[i];
    buf[j+4] = numneigh[i];
    numNbrs = nint(buf[j+4]);

    for (k = 5; k < 5+numNbrs; k++) {
      buf[j+k] = neighid[i][k-5];
    }
    j += (5+numNbrs);

    if (atom->molecule == NULL ) buf[j] = 0.0;
    else buf[j] = atom->molecule[i];
    j ++;
    
    for (k = 0; k < numNbrs; k++) {
      buf[j+k] = neigh_d[i][k];
    }
    j += (1+numNbrs);
  }
  nbuf_local = j - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::RecvBuffer(double *buf, int nbuf, int nbuf_local,
                               int natoms, int maxnum)
{
  int i, j, k, itype;
  int inode, nlocal_tmp, numNbrs;
  tagint itag,jtag;
  int nlocal = atom->nlocal;
  bigint ntimestep = update->ntimestep;
  double sbotmp, nlptmp, avqtmp, dtmp;

  MPI_Request irequest, irequest2;

  if (me == 0 ){
    fprintf(fp,"# Timestep " BIGINT_FORMAT " \n",ntimestep);
    fprintf(fp,"# \n");
    fprintf(fp,"# Number of particles %d \n",natoms);
    fprintf(fp,"# \n");
    fprintf(fp,"# Max number of bonds per atom %d\n",maxnum);
    fprintf(fp,"# Particle connection table \n");
    fprintf(fp,"# id type nb id_1...id_nb mol nlp q \n");
  }

  j = 2;
  if (me == 0) {
    for (inode = 0; inode < nprocs; inode ++) {
      if (inode == 0) {
        nlocal_tmp = nlocal;
      } else {
        MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Wait(&irequest,MPI_STATUS_IGNORE);
        nlocal_tmp = nint(buf[0]);
      }
      j = 2;
      for (i = 0; i < nlocal_tmp; i ++) {
        itag = static_cast<tagint> (buf[j-1]);
        itype = nint(buf[j+0]);
        sbotmp = buf[j+1];
        nlptmp = buf[j+2];
        avqtmp = buf[j+3];
        numNbrs = nint(buf[j+4]);

        fprintf(fp," " TAGINT_FORMAT " %d %d",itag,itype,numNbrs);

        for (k = 5; k < 5+numNbrs; k++) {
          jtag = static_cast<tagint> (buf[j+k]);
          fprintf(fp," " TAGINT_FORMAT,jtag);
        }
        j += (5+numNbrs);

        fprintf(fp," " TAGINT_FORMAT,static_cast<tagint> (buf[j]));
        j ++;
        
        for (k = 0; k < numNbrs; k++) {
          dtmp = buf[j+k];
          fprintf(fp,"%14.3f",dtmp);
        }
        j += (1+numNbrs);
        
        fprintf(fp,"%14.3f%14.3f%14.3f\n",sbotmp,nlptmp,avqtmp);
      }
    }
  } else {
    MPI_Isend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world,&irequest2);
    MPI_Wait(&irequest2,MPI_STATUS_IGNORE);
  }
  if(me ==0) fprintf(fp,"# \n");
  OfekFunc();//**********************calling my func

}

/* ---------------------------------------------------------------------- */

//TODO-ןmprove efficiency OR find builded method.


int FixReaxCOfek::from_tag_to_i(tagint tag){
  //printf("\n--------------> in from_tag_to_i function <--------------\n ") ; 
  int len=sizeof(atom->tag)/sizeof(*atom->tag);
  for(int i=0; i<len; i++){
    if(atom->tag[i]==tag){
        //printf("\nthe tag is: %d the i is: %d\n", tag, i);
        return i;
    }
      
  }
  
//printf("no i for that shitti tag");
  return -1;

}




/* ---------------------------------------------------------------------- */

void FixReaxCOfek::OfekFunc(){

  int nlocal = atom->nlocal;
printf("\n=================================\n\n");
  printf("in check bonds func");



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


     /*==========DEBUGGING==========
     printf("i-tag=%d\n", a_tag);
     printf("i=%d\n", a);
     from_tag_to_i(a_tag);
     =============================*/

      if(a_type==FIRST_TYPE){
        a_numNbr = nint(numneigh[a]);
        for (b = 0; b < a_numNbr; b++) {
          b_tag = neighid[a][b];
          i=from_tag_to_i(b_tag);
          printf("\ndistance between neigh=%f in [%d] [%d] \n", neigh_d[a][b], a, b);
          b_type = atom->type[i]; //PROBLEM!!!!! need to get b from id
              
          if(b_type==SECOND_TYPE){
              b_numNbr=nint(numneigh[i]); //PROBLEM!!!!! need to get b from id
              printf("\nfirst cond OK\n");//**********
              for (c = 0; c < b_numNbr; c++) {
                c_tag = neighid[i][c]; //PROBLEM!!!!! need to get b from id
                i=from_tag_to_i(c_tag);
                c_type = atom->type[i]; //PROBLEM!!!!! need to get type from id
                  
                    if(c_type==THIRD_TYPE){
                      printf("\nsecond cond OK\n");//**********
                      c_numNbr=nint(numneigh[i]); //PROBLEM!!!!! need to get c from id
                      for (d = 0; d < c_numNbr; d++) {
                        d_tag = neighid[i][d]; //PROBLEM!!!!! need to get c from id
                        i=from_tag_to_i(d_tag);
                        d_type = atom->type[i]; //PROBLEM!!!!! need to get type from id
                        
                          if(d_type==FORTH_TYPE){
                            printf("\nthird cond OK\n");//**********
                            d_numNbr=nint(numneigh[i]); //PROBLEM!!!!! need to get d from id
                            for (int e = 0; e < d_numNbr; e++) {
                              //PROBLEM!!!!! need to get c from id
                              if(neighid[i][e]==a_tag) {
                                
                                  printf("\n\n\n-------------->doing all the shiti math<--------------\n\n\n");
                                  /*==============
                                  ================*/
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
printf("\n\n=================================\n");
}

/* ---------------------------------------------------------------------- */

int FixReaxCOfek::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::destroy()
{
  memory->destroy(neigh_d);//**************
  memory->destroy(neighid);
  memory->destroy(numneigh);
}

/* ---------------------------------------------------------------------- */

void FixReaxCOfek::allocate()
{
  memory->create(neigh_d,nmax,MAXREAXBOND,"reax/c/ofek:neigh_d");//***************
  memory->create(neighid,nmax,MAXREAXBOND,"reax/c/ofek:neighid");
  memory->create(numneigh,nmax,"reax/c/ofek:numneigh");
}

/* ---------------------------------------------------------------------- */

double FixReaxCOfek::memory_usage()
{
  double bytes;

  bytes = 3.0*nmax*sizeof(double);
  bytes += nmax*sizeof(int);
  bytes += 1.0*nmax*MAXREAXBOND*sizeof(double);
  bytes += 1.0*nmax*MAXREAXBOND*sizeof(int);

  return bytes;
}


