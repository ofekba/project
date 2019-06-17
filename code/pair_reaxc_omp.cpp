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
   Contributing author:
   Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

   Per-atom energy/virial added by Ray Shan (Materials Design, Inc.)
   Fix reax/c/bonds and fix reax/c/species for pair_style reax/c added
   by Ray Shan (Materials Design)

   OpenMP based threading support for pair_style reax/c/omp added
   by Hasan Metin Aktulga (MSU), Chris Knight (ALCF), Paul Coffman (ALCF),
   Kurt O'Hearn (MSU), Ray Shan (Materials Design), Wei Jiang (ALCF)

   Integration of the pair_style reax/c/omp into the User-OMP package
   by Axel Kohlmeyer (Temple U.)

   Please cite the related publication:
   H. M. Aktulga, C. Knight, P. Coffman, K. A. O'Hearn, T. R. Shan,
   W. Jiang, "Optimizing the performance of reactive molecular dynamics
   simulations for multi-core architectures", International Journal of
   High Performance Computing Applications, to appear.
 ------------------------------------------------------------------------- */

#include "pair_reaxc_omp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "fix.h"
#include "fix_reaxc.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "timer.h"

#include "reaxc_types.h"
#include "reaxc_allocate.h"
#include "reaxc_control.h"
#include "reaxc_ffield.h"
#include "reaxc_forces_omp.h"
#include "reaxc_init_md_omp.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_reset_tools.h"
#include "reaxc_tool_box.h"
#include "reaxc_traj.h"
#include "reaxc_vector.h"
#include "fix_reaxc_bonds.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#ifdef OMP_TIMING
double ompTimingData[LASTTIMINGINDEX];
int ompTimingCount[LASTTIMINGINDEX];
int ompTimingCGCount[LASTTIMINGINDEX];
#endif

static const char cite_pair_reax_c_omp[] =
  "pair reax/c/omp and fix qeq/reax/omp command:\n\n"
  "@Article{Aktulga17,\n"
  " author =  {H. M. Aktulga, C. Knight, P. Coffman, K. A. OHearn, T. R. Shan, W. Jiang},\n"
  " title =   {Optimizing the performance of reactive molecular dynamics simulations for multi-core architectures},\n"
  " journal = {International Journal of High Performance Computing Applications},\n"
  " year =    to appear\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairReaxCOMP::PairReaxCOMP(LAMMPS *lmp) : PairReaxC(lmp), ThrOMP(lmp, THR_PAIR)
{
  //OFEK
  printf("\n*******in PairReaxCOMP*****\n");
  if (lmp->citeme) lmp->citeme->add(cite_pair_reax_c_omp);

  suffix_flag |= Suffix::OMP;
  system->pair_ptr = this;
  system->omp_active = 1;

  num_nbrs_offset = NULL;

#ifdef OMP_TIMING
  for (int i=0;i<LASTTIMINGINDEX;i++) {
    ompTimingData[i] = 0;
    ompTimingCount[i] = 0;
    ompTimingCGCount[i] = 0;
  }
#endif
}

/* ---------------------------------------------------------------------- */

PairReaxCOMP::~PairReaxCOMP()
{
  //OFEK
  printf("\n\n*******out PairReaxCOMP*****\n\n");
  if (setup_flag) {
    reax_list * bonds = lists+BONDS;
    for (int i=0; i<bonds->num_intrs; ++i)
      sfree(error, bonds->select.bond_list[i].bo_data.CdboReduction, "CdboReduction");
  }
  memory->destroy(num_nbrs_offset);

#ifdef OMP_TIMING
  int myrank;

  MPI_Comm_rank(mpi_data->world,&myrank);

  // Write screen output
  if (timer->has_full() && myrank == 0 && screen) {
    fprintf(screen,"\n\nWrite_Lists    took %11.3lf seconds", ompTimingData[COMPUTEWLINDEX]);

    fprintf(screen,"\n\nCompute_Forces took %11.3lf seconds:", ompTimingData[COMPUTEINDEX]);
    fprintf(screen,"\n ->Initial Forces: %11.3lf seconds", ompTimingData[COMPUTEIFINDEX]);
    fprintf(screen,"\n ->Bond Order:     %11.3lf seconds", ompTimingData[COMPUTEBOINDEX]);
    fprintf(screen,"\n ->Atom Energy:    %11.3lf seconds", ompTimingData[COMPUTEATOMENERGYINDEX]);
    fprintf(screen,"\n ->Bond:           %11.3lf seconds", ompTimingData[COMPUTEBONDSINDEX]);
    fprintf(screen,"\n ->Hydrogen bonds: %11.3lf seconds", ompTimingData[COMPUTEHBONDSINDEX]);
    fprintf(screen,"\n ->Torsion Angles: %11.3lf seconds", ompTimingData[COMPUTETORSIONANGLESBOINDEX]);
    fprintf(screen,"\n ->Valence Angles: %11.3lf seconds", ompTimingData[COMPUTEVALENCEANGLESBOINDEX]);
    fprintf(screen,"\n ->Non-Bonded For: %11.3lf seconds", ompTimingData[COMPUTENBFINDEX]);
    fprintf(screen,"\n ->Total Forces:   %11.3lf seconds", ompTimingData[COMPUTETFINDEX]);

    fprintf(screen,"\n\nfixQEQ:          %11.3lf seconds", ompTimingData[COMPUTEQEQINDEX]);
    fprintf(screen,"\n ->QEQ init:       %11.3lf seconds", ompTimingData[COMPUTEINITMVINDEX]);

    double avg = double(ompTimingCGCount[COMPUTECG1INDEX]) / double(ompTimingCount[COMPUTECG1INDEX]);
    fprintf(screen,"\n ->QEQ CG1:        %11.3lf seconds with %4.1lf iterations on average.", ompTimingData[COMPUTECG1INDEX], avg);

    avg = double(ompTimingCGCount[COMPUTECG2INDEX]) / double(ompTimingCount[COMPUTECG2INDEX]);
    fprintf(screen,"\n ->QEQ CG2:        %11.3lf seconds with %4.1lf iterations on average.", ompTimingData[COMPUTECG2INDEX], avg);
    fprintf(screen,"\n ->QEQ CalcQ:      %11.3lf seconds\n", ompTimingData[COMPUTECALCQINDEX]);
  }

  // Write logfile output
  if (timer->has_full() && myrank == 0 && logfile) {
    fprintf(logfile,"\n\nWrite_Lists    took %11.3lf seconds", ompTimingData[COMPUTEWLINDEX]);

    fprintf(logfile,"\n\nCompute_Forces took %11.3lf seconds:", ompTimingData[COMPUTEINDEX]);
    fprintf(logfile,"\n ->Initial Forces: %11.3lf seconds", ompTimingData[COMPUTEIFINDEX]);
    fprintf(logfile,"\n ->Bond Order:     %11.3lf seconds", ompTimingData[COMPUTEBOINDEX]);
    fprintf(logfile,"\n ->Atom Energy:    %11.3lf seconds", ompTimingData[COMPUTEATOMENERGYINDEX]);
    fprintf(logfile,"\n ->Bond:           %11.3lf seconds", ompTimingData[COMPUTEBONDSINDEX]);
    fprintf(logfile,"\n ->Hydrogen bonds: %11.3lf seconds", ompTimingData[COMPUTEHBONDSINDEX]);
    fprintf(logfile,"\n ->Torsion Angles: %11.3lf seconds", ompTimingData[COMPUTETORSIONANGLESBOINDEX]);
    fprintf(logfile,"\n ->Valence Angles: %11.3lf seconds", ompTimingData[COMPUTEVALENCEANGLESBOINDEX]);
    fprintf(logfile,"\n ->Non-Bonded For: %11.3lf seconds", ompTimingData[COMPUTENBFINDEX]);
    fprintf(logfile,"\n ->Total Forces:   %11.3lf seconds", ompTimingData[COMPUTETFINDEX]);

    fprintf(logfile,"\n\nfixQEQ:          %11.3lf seconds", ompTimingData[COMPUTEQEQINDEX]);
    fprintf(logfile,"\n ->QEQ init:       %11.3lf seconds", ompTimingData[COMPUTEINITMVINDEX]);

    double avg = double(ompTimingCGCount[COMPUTECG1INDEX]) / double(ompTimingCount[COMPUTECG1INDEX]);
    fprintf(logfile,"\n ->QEQ CG1:        %11.3lf seconds with %4.1lf iterations on average.", ompTimingData[COMPUTECG1INDEX], avg);

    avg = double(ompTimingCGCount[COMPUTECG2INDEX]) / double(ompTimingCount[COMPUTECG2INDEX]);
    fprintf(logfile,"\n ->QEQ CG2:        %11.3lf seconds with %4.1lf iterations on average.", ompTimingData[COMPUTECG2INDEX], avg);
    fprintf(logfile,"\n ->QEQ CalcQ:      %11.3lf seconds\n", ompTimingData[COMPUTECALCQINDEX]);
  }
#endif
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::compute(int eflag, int vflag)
{
  double evdwl,ecoul;
  double t_start, t_end;

  //ofek
  //set_extra_potential_parameters();

  // communicate num_bonds once every reneighboring
  // 2 num arrays stored by fix, grab ptr to them

  if (neighbor->ago == 0) comm->forward_comm_fix(fix_reax);
  int *num_bonds = fix_reax->num_bonds;
  int *num_hbonds = fix_reax->num_hbonds;

  if(calm_down>0){
    calm_down--;
  }

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  if (vflag_global) control->virial = 1;
  else control->virial = 0;

  if (vflag_atom)
     error->all(FLERR,"Pair style reax/c/omp does not support "
                "computing per-atom stress");

  system->n = atom->nlocal; // my atoms
  system->N = atom->nlocal + atom->nghost; // mine + ghosts
  system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system

  system->big_box.V = 0;
  system->big_box.box_norms[0] = 0;
  system->big_box.box_norms[1] = 0;
  system->big_box.box_norms[2] = 0;
  if (comm->me == 0 ) t_start = MPI_Wtime();
  // setup data structures

  setup();

  Reset( system, control, data, workspace, &lists );

  // Why not update workspace like in MPI-only code?
  // Using the MPI-only way messes up the hb energy
  //workspace->realloc.num_far = write_reax_lists();
  write_reax_lists();

  // timing for filling in the reax lists
  if (comm->me == 0) {
    t_end = MPI_Wtime();
    data->timing.nbrs = t_end - t_start;
  }

  // forces

#ifdef OMP_TIMING
  double startTimeBase,endTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  Compute_ForcesOMP(system,control,data,workspace,&lists,out_control,mpi_data);
  
  // mine
  double added_e=0; //the energy of the BB potential
	if(flag_bb==1){
    added_e=compute_BB();
    if(added_e==-1){
      printf("\n\tERROR! FAILED OPERATE THE EXTRA POTENTIAL\n");
      added_e=0;
    }
    eng_vdwl += added_e;
  }
  if(energy_fp!=NULL)
    fprintf(energy_fp,"\n%f",added_e);
  else
    printf("cant write\n");
    //OFEK ^^

  read_reax_forces(vflag);

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEINDEX] += (endTimeBase-startTimeBase);
#endif

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for(int k = 0; k < system->N; ++k) {
    num_bonds[k] = system->my_atoms[k].num_bonds;
    num_hbonds[k] = system->my_atoms[k].num_hbonds;
  }

  // energies and pressure

  if (eflag_global) {
    evdwl += data->my_en.e_bond;
    evdwl += data->my_en.e_ov;
    evdwl += data->my_en.e_un;
    evdwl += data->my_en.e_lp;
    evdwl += data->my_en.e_ang;
    evdwl += data->my_en.e_pen;
    evdwl += data->my_en.e_coa;
    evdwl += data->my_en.e_hb;
    evdwl += data->my_en.e_tor;
    evdwl += data->my_en.e_con;
    evdwl += data->my_en.e_vdW;

    ecoul += data->my_en.e_ele;
    ecoul += data->my_en.e_pol;

    // Store the different parts of the energy
    // in a list for output by compute pair command

    pvector[0] = data->my_en.e_bond;
    pvector[1] = data->my_en.e_ov + data->my_en.e_un;
    pvector[2] = data->my_en.e_lp;
    pvector[3] = 0.0;
    pvector[4] = data->my_en.e_ang;
    pvector[5] = data->my_en.e_pen;
    pvector[6] = data->my_en.e_coa;
    pvector[7] = data->my_en.e_hb;
    pvector[8] = data->my_en.e_tor;
    pvector[9] = data->my_en.e_con;
    pvector[10] = data->my_en.e_vdW;
    pvector[11] = data->my_en.e_ele;
    pvector[12] = 0.0;
    pvector[13] = data->my_en.e_pol;
  }

  if (vflag_fdotr) virial_fdotr_compute();

// Set internal timestep counter to that of LAMMPS

  data->step = update->ntimestep;

  Output_Results( system, control, data, &lists, out_control, mpi_data );

  // populate tmpid and tmpbo arrays for fix reax/c/species
  int i, j;

  if(fixspecies_flag) {
    if (system->N > nmax) {
      memory->destroy(tmpid);
      memory->destroy(tmpbo);
      nmax = system->N;
      memory->create(tmpid,nmax,MAXSPECBOND,"pair:tmpid");
      memory->create(tmpbo,nmax,MAXSPECBOND,"pair:tmpbo");
    }

#if defined(_OPENMP)
#pragma omp parallel for collapse(2) schedule(static) default(shared)
#endif
    for (i = 0; i < system->N; i ++)
      for (j = 0; j < MAXSPECBOND; j ++) {
        tmpbo[i][j] = 0.0;
        tmpid[i][j] = 0;
      }

    FindBond();
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::init_style( )
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair reax/c/omp requires atom attribute q");

  // firstwarn = 1;

  int iqeq = modify->find_fix_by_style("qeq/reax/omp");
  if (iqeq < 0 && qeqflag == 1)
    error->all(FLERR,"Pair reax/c/omp requires use of fix qeq/reax/omp");

  system->n = atom->nlocal; // my atoms
  system->N = atom->nlocal + atom->nghost; // mine + ghosts
  system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system
  system->wsize = comm->nprocs;

  system->big_box.V = 0;
  system->big_box.box_norms[0] = 0;
  system->big_box.box_norms[1] = 0;
  system->big_box.box_norms[2] = 0;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style reax/c/omp requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style reax/c/omp requires newton pair on");

 if ((atom->map_tag_max > 99999999) && (comm->me == 0))
    error->warning(FLERR,"Some Atom-IDs are too large. Pair style reax/c/omp "
                   "native output files may get misformatted or corrupted");

  // because system->bigN is an int, we cannot have more atoms than MAXSMALLINT

  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Too many atoms for pair style reax/c/omp");

  // need a half neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;

  cutmax = MAX3(control->nonb_cut, control->hbond_cut, control->bond_cut);
  if ((cutmax < 2.0*control->bond_cut) && (comm->me == 0))
    error->warning(FLERR,"Total cutoff < 2*bond cutoff. May need to use an "
                   "increased neighbor list skin.");

  for( int i = 0; i < LIST_N; ++i )
    lists[i].allocated = 0;

  if (fix_reax == NULL) {
    char **fixarg = new char*[3];
    fixarg[0] = (char *) fix_id;
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "REAXC";
    modify->add_fix(3,fixarg);
    delete [] fixarg;
    fix_reax = (FixReaxC *) modify->fix[modify->nfix-1];
  }

#if defined(_OPENMP)
  control->nthreads = omp_get_max_threads();
#else
  control->nthreads = 1;
#endif
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::setup( )
{
  int oldN;
  int mincap = system->mincap;
  double safezone = system->safezone;

  system->n = atom->nlocal; // my atoms
  system->N = atom->nlocal + atom->nghost; // mine + ghosts
  oldN = system->N;
  system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system

  if (system->N > nmax) {
    memory->destroy(num_nbrs_offset);
    // Don't update nmax here. It is updated at end of compute().
    memory->create(num_nbrs_offset, system->N, "pair:num_nbrs_offset");
  }

  if (setup_flag == 0) {

    setup_flag = 1;

    int *num_bonds = fix_reax->num_bonds;
    int *num_hbonds = fix_reax->num_hbonds;

    control->vlist_cut = neighbor->cutneighmax;

    // determine the local and total capacity

    system->local_cap = MAX( (int)(system->n * safezone), mincap );
    system->total_cap = MAX( (int)(system->N * safezone), mincap );

    // initialize my data structures

    PreAllocate_Space( system, control, workspace );
    write_reax_atoms();

    int num_nbrs = estimate_reax_lists();
    if(!Make_List(system->total_cap, num_nbrs, TYP_FAR_NEIGHBOR,
                  lists+FAR_NBRS))
      error->all(FLERR,"Pair reax/c problem in far neighbor list");

    write_reax_lists();

    InitializeOMP( system, control, data, workspace, &lists, out_control,
                mpi_data, world );

    for( int k = 0; k < system->N; ++k ) {
      num_bonds[k] = system->my_atoms[k].num_bonds;
      num_hbonds[k] = system->my_atoms[k].num_hbonds;
    }

  } else {

    // fill in reax datastructures

    write_reax_atoms();

    // reset the bond list info for new atoms

    for(int k = oldN; k < system->N; ++k)
      Set_End_Index( k, Start_Index( k, lists+BONDS ), lists+BONDS );

    // estimate far neighbor list size
    // Not present in MPI-only version
    workspace->realloc.num_far = estimate_reax_lists();

    // check if I need to shrink/extend my data-structs

    ReAllocate( system, control, data, workspace, &lists );
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::write_reax_atoms()
{
  int *num_bonds = fix_reax->num_bonds;
  int *num_hbonds = fix_reax->num_hbonds;

  if (system->N > system->total_cap)
    error->all(FLERR,"Too many ghost atoms");

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)
#endif
  for( int i = 0; i < system->N; ++i ){
    system->my_atoms[i].orig_id = atom->tag[i];
    system->my_atoms[i].type = map[atom->type[i]];
    system->my_atoms[i].x[0] = atom->x[i][0];
    system->my_atoms[i].x[1] = atom->x[i][1];
    system->my_atoms[i].x[2] = atom->x[i][2];
    system->my_atoms[i].q = atom->q[i];
    system->my_atoms[i].num_bonds = num_bonds[i];
    system->my_atoms[i].num_hbonds = num_hbonds[i];
  }
}

/* ---------------------------------------------------------------------- */

int PairReaxCOMP::estimate_reax_lists()
{
  int i;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int numall = list->inum + list->gnum;
  int mincap = system->mincap;

  // for good performance in the OpenMP implementation, each thread needs
  // to know where to place the neighbors of the atoms it is responsible for.
  // The sumscan values for the list->numneigh will be used to determine the
  // neighbor offset of each atom. Note that this may cause some significant
  // memory overhead if delayed neighboring is used - so it may be desirable
  // to work on this part to reduce the memory footprint of the far_nbrs list.

  int num_nbrs = 0;

  for (int itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    num_nbrs += numneigh[i];
  }

  int new_estimate = MAX (num_nbrs, mincap*MIN_NBRS);

  return new_estimate;
}

/* ---------------------------------------------------------------------- */

int PairReaxCOMP::write_reax_lists()
{
#ifdef OMP_TIMING
  double startTimeBase, endTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  int itr_i, itr_j, i, j, num_mynbrs;
  int *jlist;
  double d_sqr, dist, cutoff_sqr;
  rvec dvec;

  double **x = atom->x;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  reax_list *far_nbrs = lists + FAR_NBRS;
  far_neighbor_data *far_list = far_nbrs->select.far_nbr_list;

  int num_nbrs = 0;
  int inum = list->inum;
  int gnum = list->gnum;
  int numall = inum + gnum;

  // sumscan of the number of neighbors per atom to determine the offsets
  // most likely, we are overallocating. desirable to work on this part
  // to reduce the memory footprint of the far_nbrs list.

  num_nbrs = 0;

  for (itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    num_nbrs_offset[i] = num_nbrs;
    num_nbrs += numneigh[i];
  }

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) default(shared)           \
  private(itr_i, itr_j, i, j, jlist, cutoff_sqr, num_mynbrs, d_sqr, dvec, dist)
#endif
  for (itr_i = 0; itr_i < numall; ++itr_i) {
    i = ilist[itr_i];
    jlist = firstneigh[i];
    Set_Start_Index( i, num_nbrs_offset[i], far_nbrs );

    if (i < inum)
      cutoff_sqr = control->nonb_cut*control->nonb_cut;
    else
      cutoff_sqr = control->bond_cut*control->bond_cut;

    num_mynbrs = 0;

    for (itr_j = 0; itr_j < numneigh[i]; ++itr_j) {
      j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance( x[j], x[i], &d_sqr, &dvec );

      if (d_sqr <= cutoff_sqr) {
        dist = sqrt( d_sqr );
        set_far_nbr( &far_list[num_nbrs_offset[i] + num_mynbrs], j, dist, dvec );
        ++num_mynbrs;
      }
    }
    Set_End_Index( i, num_nbrs_offset[i] + num_mynbrs, far_nbrs );
  }

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEWLINDEX] += (endTimeBase-startTimeBase);
#endif

  return num_nbrs;
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::read_reax_forces(int /* vflag */)
{
  if(flag_bb==1){
    if(count_bb_timesteps<MAX_NUM_TIMESTEPS){
      count_bb_timesteps++;
    }
    else{
      flag_bb=0;
      count_bb_timesteps=0;
      if(energy_fp!=NULL)
        fprintf(energy_fp,"\nfinish");
      calm_down=1000;
      printf("\n\n**** finish %d timesteps at timestep %d****\n\n", MAX_NUM_TIMESTEPS, update->ntimestep); 
    }
  }


#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)
#endif
  for( int i = 0; i < system->N; ++i ) {
    system->my_atoms[i].f[0] = workspace->f[i][0];
    system->my_atoms[i].f[1] = workspace->f[i][1];
    system->my_atoms[i].f[2] = workspace->f[i][2];

    atom->f[i][0] = -workspace->f[i][0];
    atom->f[i][1] = -workspace->f[i][1];
    atom->f[i][2] = -workspace->f[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void PairReaxCOMP::FindBond()
{
  const double bo_cut = 0.10;
  int i;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) default(shared)   \
  private(i)
#endif
  for (i = 0; i < system->n; i++) {
    int j, pj, nj;
    double bo_tmp;
    bond_data *bo_ij;

    nj = 0;
    for( pj = Start_Index(i, lists); pj < End_Index(i, lists); ++pj ) {
      bo_ij = &( lists->select.bond_list[pj] );
      j = bo_ij->nbr;
      if (j < i) continue;

      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp >= bo_cut ) {
        tmpid[i][nj] = j;
        tmpbo[i][nj] = bo_tmp;
        nj ++;
        if (nj > MAXSPECBOND) error->all(FLERR,"Increase MAXSPECBOND in fix_reaxc_species.h");
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
/*
double PairReaxCOMP::compute_BB(){
  //printf("\nin computebb ()\n");
  for (int i = 0; i < atom->nlocal; i++) {
    for (int j = 0; j < 3; j++) {
      f_fourset[i][j] = 0;
    }
  }
  double e=0; //the amount of the additional energy

  for(int k=0; k<num_fourset; k++){
    e+=compute_BB_pair(fourset[k][0], fourset[k][1]); //O-H
    if(e==-1) return-1;
    e+=compute_BB_pair(fourset[k][0], fourset[k][3]); //O-C
    if(e==-1) return-1;
    e+=compute_BB_pair(fourset[k][2], fourset[k][3]); //N-C
    if(e==-1) return-1;
    //e+=compute_BB_pair(fourset[k][2], fourset[k][1]); //N-H
    //if(e==-1) return-1;
  }

  //printf("\nout computebb ()\n");
 return e;
} */

/* ---------------------------------------------------------------------- */

/* double PairReaxCOMP::compute_BB_pair(int i_tag, int j_tag){
 // if(update->ntimestep>1000 || update->ntimestep<1010) printf("\n\n\nin PairReaxCOMP::compute_BB_pair\n\n\n");
  
  int i,j,itype,jtype,pk,k;
  double fpair, rij=0;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int cal_dist_flag=0;
  double d_sqr;
  rvec dvec, xi, xj;
  
  i=tag_to_i(i_tag);
  if(i==-1) return -1;
  if(atom->tag[i]!=i_tag) return -1;
  itype = type[i];
  
  j=tag_to_i(j_tag);
  if(j==-1) return -1;
  if(atom->tag[j]!=j_tag) return -1;
  jtype = type[j];

  bond_data *bo_ij;
  #if defined(_OPENMP)
  #pragma omp parallel for schedule(static) default(shared)
  #endif
  for( pk = Start_Index(i, lists+BONDS); pk < End_Index(i, lists+BONDS); ++pk ) {
      bo_ij = &( (lists+BONDS)->select.bond_list[pk] );
      k = bo_ij->nbr;
      if (j == k || atom->tag[k] == j_tag){
        cal_dist_flag=1;
        rij=bo_ij->d;
        dvec[0]=bo_ij->dvec[0];
        dvec[1]=bo_ij->dvec[1];
        dvec[2]=bo_ij->dvec[2];
       // break;
      }
  }

  if(rij==0){
    far_neighbor_data *nbr_ij;
  #if defined(_OPENMP)
  #pragma omp parallel for schedule(static) default(shared)
  #endif
    for( pk = Start_Index(i, lists+FAR_NBRS); pk < End_Index(i, lists+FAR_NBRS); ++pk ) {
      nbr_ij = &( (lists+FAR_NBRS)->select.far_nbr_list[pk] );
      k = nbr_ij->nbr;
      if (j == k || atom->tag[k] == j_tag){
        cal_dist_flag=2;
        rij=nbr_ij->d;
        dvec[0]=nbr_ij->dvec[0];
        dvec[1]=nbr_ij->dvec[1];
        dvec[2]=nbr_ij->dvec[2];
       // break;
      }
    }
  }
  
  if(rij==0){
    xi[0]=x[i][0];
    xi[1]=x[i][1];
    xi[2]=x[i][2];
    xj[0]=x[j][0];
    xj[1]=x[j][1];
    xj[2]=x[j][2];
    get_distance(xj, xi, &d_sqr, &dvec );
    rij=sqrt(d_sqr);
    cal_dist_flag=3;
  }
  int sign;
  if((x[i][0]-x[j][0]) * dvec[0] > 0)  sign=-1;
  else  sign=1;
    
  //OFEK
  //FOR DEBUGGING
  //print the distance at the first 2 timesteps and at the last 2 timesteps
  //if(count_bb_timesteps<1 || count_bb_timesteps>MAX_NUM_TIMESTEPS-1 ){
  if(count_bb_timesteps%1000==0 ){
  	printf("\ncal_dist_flag=%d",cal_dist_flag);
    if( (itype==1 && jtype==4))
      printf("\nThe distance between C (TAG=%d) ,N(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==4 && jtype==1))
      printf("\nThe distance between N (TAG=%d) ,C(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==3 && jtype==2))
      printf("\nThe distance between O (TAG=%d) ,H(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==2 && jtype==3))
      printf("\nThe distance between H (TAG=%d) ,O(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==3 && jtype==1))
      printf("\nThe distance between O (TAG=%d) ,C(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==1 && jtype==3))
      printf("\nThe distance between C (TAG=%d) ,O(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==2 && jtype==4))
      printf("\nThe distance between H (TAG=%d) ,N(TAG=%d) =%f\n", i_tag, j_tag, rij);
    else if( (itype==4 && jtype==2))
      printf("\nThe distance between N (TAG=%d) ,H(TAG=%d) =%f\n", i_tag, j_tag, rij);
  }
    
  fpair=single_BB(i, j, i_tag, j_tag, itype, jtype, rij);
    
  // //calculate the F force vector for atom i
  // f_fourset[i_tag-1][0] += (delx*fpair)/rij;
  // f_fourset[i_tag-1][1] += (dely*fpair)/rij;
  // f_fourset[i_tag-1][2] += (delz*fpair)/rij;
  // //calculate the F force vector for atom j       
  // f_fourset[j_tag-1][0] -= (delx*fpair)/rij;
  // f_fourset[j_tag-1][1] -= (dely*fpair)/rij;
  // f_fourset[j_tag-1][2] -= (delz*fpair)/rij;

  workspace->f[i][0] -= sign*(dvec[0]*fpair)/rij;
  workspace->f[i][1] -= sign*(dvec[1]*fpair)/rij;
  workspace->f[i][2] -= sign*(dvec[2]*fpair)/rij;

  workspace->f[j][0] += sign*(dvec[0]*fpair)/rij;
  workspace->f[j][1] += sign*(dvec[1]*fpair)/rij;
  workspace->f[j][2] += sign*(dvec[2]*fpair)/rij;
    
  //calculate and return the E (energy)
  double r=rij-wanted_dist[itype][jtype];
  double e=F1[itype][jtype] * (1 - exp( -F2[itype][jtype] * r * r ));
  return e;

}*/