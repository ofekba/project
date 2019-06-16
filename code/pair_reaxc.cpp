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
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)
   Per-atom energy/virial added by Ray Shan (Sandia)
   Fix reax/c/bonds and fix reax/c/species for pair_style reax/c added by
        Ray Shan (Sandia)
   Hybrid and hybrid/overlay compatibility added by Ray Shan (Sandia)
------------------------------------------------------------------------- */

#include "pair_reaxc.h"
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

#include "reaxc_types.h"
#include "reaxc_allocate.h"
#include "reaxc_control.h"
#include "reaxc_ffield.h"
#include "reaxc_forces.h"
#include "reaxc_init_md.h"
#include "reaxc_io_tools.h"
#include "reaxc_list.h"
#include "reaxc_lookup.h"
#include "reaxc_reset_tools.h"
#include "reaxc_traj.h"
#include "reaxc_vector.h"
#include "fix_reaxc_bonds.h"

using namespace LAMMPS_NS;

static const char cite_pair_reax_c[] =
  "pair reax/c command:\n\n"
  "@Article{Aktulga12,\n"
  " author = {H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama},\n"
  " title = {Parallel reactive molecular dynamics: Numerical methods and algorithmic techniques},\n"
  " journal = {Parallel Computing},\n"
  " year =    2012,\n"
  " volume =  38,\n"
  " pages =   {245--259}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairReaxC::PairReaxC(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_reax_c);

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  ghostneigh = 1;

  fix_id = new char[24];
  snprintf(fix_id,24,"REAXC_%d",instance_me);

  system = (reax_system *)
    memory->smalloc(sizeof(reax_system),"reax:system");
  memset(system,0,sizeof(reax_system));
  control = (control_params *)
    memory->smalloc(sizeof(control_params),"reax:control");
  memset(control,0,sizeof(control_params));
  data = (simulation_data *)
    memory->smalloc(sizeof(simulation_data),"reax:data");
  workspace = (storage *)
    memory->smalloc(sizeof(storage),"reax:storage");
  lists = (reax_list *)
    memory->smalloc(LIST_N * sizeof(reax_list),"reax:lists");
  memset(lists,0,LIST_N * sizeof(reax_list));
  out_control = (output_controls *)
    memory->smalloc(sizeof(output_controls),"reax:out_control");
  memset(out_control,0,sizeof(output_controls));
  mpi_data = (mpi_datatypes *)
    memory->smalloc(sizeof(mpi_datatypes),"reax:mpi");

  control->me = system->my_rank = comm->me;

  system->my_coords[0] = 0;
  system->my_coords[1] = 0;
  system->my_coords[2] = 0;
  system->num_nbrs = 0;
  system->n = 0; // my atoms
  system->N = 0; // mine + ghosts
  system->bigN = 0;  // all atoms in the system
  system->local_cap = 0;
  system->total_cap = 0;
  system->gcell_cap = 0;
  system->bndry_cuts.ghost_nonb = 0;
  system->bndry_cuts.ghost_hbond = 0;
  system->bndry_cuts.ghost_bond = 0;
  system->bndry_cuts.ghost_cutoff = 0;
  system->my_atoms = NULL;
  system->pair_ptr = this;
  system->error_ptr = error;
  control->error_ptr = error;

  system->omp_active = 0;

  fix_reax = NULL;
  tmpid = NULL;
  tmpbo = NULL;

  nextra = 14;
  pvector = new double[nextra];

  setup_flag = 0;
  fixspecies_flag = 0;

  nmax = 0;
  
  //mine
  f_fourset=NULL; //ofek no need
  fourset=NULL;
  num_fourset=0;
  count_bb_timesteps=0;
  flag_bb=0;
  wanted_dist=NULL;
  F1=NULL;
  F2=NULL;
  MAX_NUM_TIMESTEPS=5000; //DEFAULT VALUE
  parameters_fp=NULL;
  calm_down=0;


  //FOR ENERGY_FP
  energy_fp = fopen("energy.reax","w");
  if (energy_fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open fix reax/c/bonds file energy.reax");
    error->one(FLERR,str);
  }



}

/* ---------------------------------------------------------------------- */

PairReaxC::~PairReaxC()
{
  printf("\n\n*******out PairReaxC*****\n\n");
  if (copymode) return;

  if (fix_reax) modify->delete_fix(fix_id);
  delete[] fix_id;

  if (setup_flag) {
    Close_Output_Files( system, control, out_control, mpi_data );

    // deallocate reax data-structures

    if (control->tabulate ) Deallocate_Lookup_Tables( system);

    if (control->hbond_cut > 0 )  Delete_List( lists+HBONDS );
    Delete_List( lists+BONDS );
    Delete_List( lists+THREE_BODIES );
    Delete_List( lists+FAR_NBRS );

    DeAllocate_Workspace( control, workspace );
    DeAllocate_System( system );
  }

  memory->destroy( system );
  memory->destroy( control );
  memory->destroy( data );
  memory->destroy( workspace );
  memory->destroy( lists );
  memory->destroy( out_control );
  memory->destroy( mpi_data );

  //mine
  memory->destroy( f_fourset );
  memory->destroy( fourset );
  memory->destroy( wanted_dist );
  memory->destroy( F1 );
  memory->destroy( F2 );
  fclose(energy_fp);

  // deallocate interface storage
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    delete [] map;

    delete [] chi;
    delete [] eta;
    delete [] gamma;
  }

  memory->destroy(tmpid);
  memory->destroy(tmpbo);

  delete [] pvector;

}

/* ---------------------------------------------------------------------- */

void PairReaxC::allocate( )
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  //mine
  memory->create(f_fourset,atom->nlocal,3,"pair:f_fourset"); //no need
  memory->create(fourset,atom->nlocal,4,"pair:fourset");
  memory->create(wanted_dist,n+1,n+1,"pair:wanted_dist");
  memory->create(F1,n+1,n+1,"pair:F1");
  memory->create(F2,n+1,n+1,"pair:F2");

  map = new int[n+1];

  chi = new double[n+1];
  eta = new double[n+1];
  gamma = new double[n+1];
}

/* ---------------------------------------------------------------------- */

void PairReaxC::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");

  // read name of control file or use default controls

  if (strcmp(arg[0],"NULL") == 0) {
    strcpy( control->sim_name, "simulate" );
    control->ensemble = 0;
    out_control->energy_update_freq = 0;
    control->tabulate = 0;

    control->reneighbor = 1;
    control->vlist_cut = control->nonb_cut;
    control->bond_cut = 5.;
    control->hbond_cut = 7.50;
    control->thb_cut = 0.001;
    control->thb_cutsq = 0.00001;
    control->bg_cut = 0.3;

    // Initialize for when omp style included
    control->nthreads = 1;

    out_control->write_steps = 0;
    out_control->traj_method = 0;
    strcpy( out_control->traj_title, "default_title" );
    out_control->atom_info = 0;
    out_control->bond_info = 0;
    out_control->angle_info = 0;
  } else Read_Control_File(arg[0], control, out_control);

  // default values

  qeqflag = 1;
  control->lgflag = 0;
  control->enobondsflag = 1;
  system->mincap = MIN_CAP;
  system->safezone = SAFE_ZONE;
  system->saferzone = SAFER_ZONE;

  // process optional keywords

  int iarg = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"checkqeq") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
      if (strcmp(arg[iarg+1],"yes") == 0) qeqflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) qeqflag = 0;
      else error->all(FLERR,"Illegal pair_style reax/c command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"enobonds") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
      if (strcmp(arg[iarg+1],"yes") == 0) control->enobondsflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) control->enobondsflag = 0;
      else error->all(FLERR,"Illegal pair_style reax/c command");
      iarg += 2;
  } else if (strcmp(arg[iarg],"lgvdw") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
      if (strcmp(arg[iarg+1],"yes") == 0) control->lgflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) control->lgflag = 0;
      else error->all(FLERR,"Illegal pair_style reax/c command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"safezone") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
      system->safezone = force->numeric(FLERR,arg[iarg+1]);
      if (system->safezone < 0.0)
        error->all(FLERR,"Illegal pair_style reax/c safezone command");
      system->saferzone = system->safezone*1.2 + 0.2;
      iarg += 2;
    } else if (strcmp(arg[iarg],"mincap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style reax/c command");
      system->mincap = force->inumeric(FLERR,arg[iarg+1]);
      if (system->mincap < 0)
        error->all(FLERR,"Illegal pair_style reax/c mincap command");
      iarg += 2;
    } else error->all(FLERR,"Illegal pair_style reax/c command");
  }

  // LAMMPS is responsible for generating nbrs

  control->reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

void PairReaxC::coeff( int nargs, char **args )
{
  if (!allocated) allocate();

  if (nargs != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(args[0],"*") != 0 || strcmp(args[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read ffield file

  char *file = args[2];
  FILE *fp;
  fp = force->open_potential(file);
  if (fp != NULL)
    Read_Force_Field(fp, &(system->reax_param), control);
  else {
      char str[128];
      snprintf(str,128,"Cannot open ReaxFF potential file %s",file);
      error->all(FLERR,str);
  }

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  int itmp = 0;
  int nreax_types = system->reax_param.num_atom_types;
  for (int i = 3; i < nargs; i++) {
    if (strcmp(args[i],"NULL") == 0) {
      map[i-2] = -1;
      itmp ++;
      continue;
    }
  }

  int n = atom->ntypes;

  // pair_coeff element map
  for (int i = 3; i < nargs; i++)
    for (int j = 0; j < nreax_types; j++)
      if (strcasecmp(args[i],system->reax_param.sbp[j].name) == 0) {
        map[i-2] = j;
        itmp ++;
      }

  // error check
  if (itmp != n)
    error->all(FLERR,"Non-existent ReaxFF type");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  //mine
  //not here
  for(int i=0; i<atom->ntypes+1; i++){
    for(int j=0; j<atom->ntypes+1; j++){
      F1[i][j]=F2[i][j]=0;
      wanted_dist[i][j]=0;
    }
  }

}

/* ---------------------------------------------------------------------- */

void PairReaxC::init_style( )
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style reax/c requires atom attribute q");

  // firstwarn = 1;

  int iqeq;
  for (iqeq = 0; iqeq < modify->nfix; iqeq++)
    if (strstr(modify->fix[iqeq]->style,"qeq/reax")
       || strstr(modify->fix[iqeq]->style,"qeq/shielded")) break;
  if (iqeq == modify->nfix && qeqflag == 1)
    error->all(FLERR,"Pair reax/c requires use of fix qeq/reax");

  system->n = atom->nlocal; // my atoms
  system->N = atom->nlocal + atom->nghost; // mine + ghosts
  system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system
  system->wsize = comm->nprocs;

  system->big_box.V = 0;
  system->big_box.box_norms[0] = 0;
  system->big_box.box_norms[1] = 0;
  system->big_box.box_norms[2] = 0;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style reax/c requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style reax/c requires newton pair on");
  if ((atom->map_tag_max > 99999999) && (comm->me == 0))
    error->warning(FLERR,"Some Atom-IDs are too large. Pair style reax/c "
                   "native output files may get misformatted or corrupted");

  // because system->bigN is an int, we cannot have more atoms than MAXSMALLINT

  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Too many atoms for pair style reax/c");

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
    if (lists[i].allocated != 1)
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
}

/* ---------------------------------------------------------------------- */

void PairReaxC::setup( )
{
  int oldN;
  int mincap = system->mincap;
  double safezone = system->safezone;

  system->n = atom->nlocal; // my atoms
  system->N = atom->nlocal + atom->nghost; // mine + ghosts
  oldN = system->N;
  system->bigN = static_cast<int> (atom->natoms);  // all atoms in the system

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
      error->one(FLERR,"Pair reax/c problem in far neighbor list");
    (lists+FAR_NBRS)->error_ptr=error;

    write_reax_lists();
    Initialize( system, control, data, workspace, &lists, out_control,
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

    // check if I need to shrink/extend my data-structs

    ReAllocate( system, control, data, workspace, &lists );
  }

  bigint local_ngroup = list->inum;
  MPI_Allreduce( &local_ngroup, &ngroup, 1, MPI_LMP_BIGINT, MPI_SUM, world );
}

/* ---------------------------------------------------------------------- */

double PairReaxC::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cutghost[i][j] = cutghost[j][i] = cutmax;
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairReaxC::compute(int eflag, int vflag)
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


  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  if (vflag_global) control->virial = 1;
  else control->virial = 0;

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
  workspace->realloc.num_far = write_reax_lists();
  // timing for filling in the reax lists
  if (comm->me == 0) {
    t_end = MPI_Wtime();
    data->timing.nbrs = t_end - t_start;
  }

  // forces

  Compute_Forces(system,control,data,workspace,&lists,out_control,mpi_data);

  // mine
  double added_e=0; //the energy of the BB potential
	if(flag_bb==1){
    added_e=compute_BB();
    if(added_e==-1){
      //OFEK
      printf("\n\tERROR! FAILED OPERATE THE EXTRA POTENTIAL\n");
      added_e=0;
    }
    eng_vdwl += added_e;
  }
  if(energy_fp!=NULL)
    fprintf(energy_fp,"\n%f",added_e);
  else
  printf("cant write\n");
  //OFEK^^

  read_reax_forces(vflag);

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

    // eng_vdwl += evdwl;
    // eng_coul += ecoul;

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

    for (i = 0; i < system->N; i ++)
      for (j = 0; j < MAXSPECBOND; j ++) {
        tmpbo[i][j] = 0.0;
        tmpid[i][j] = 0;
      }
    FindBond();
  }

}

/* ---------------------------------------------------------------------- */


void PairReaxC::write_reax_atoms()
{
  int *num_bonds = fix_reax->num_bonds;
  int *num_hbonds = fix_reax->num_hbonds;

  if (system->N > system->total_cap)
    error->all(FLERR,"Too many ghost atoms");

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

void PairReaxC::get_distance( rvec xj, rvec xi, double *d_sqr, rvec *dvec )
{
  (*dvec)[0] = xj[0] - xi[0];
  (*dvec)[1] = xj[1] - xi[1];
  (*dvec)[2] = xj[2] - xi[2];
  *d_sqr = SQR((*dvec)[0]) + SQR((*dvec)[1]) + SQR((*dvec)[2]);
}

/* ---------------------------------------------------------------------- */

void PairReaxC::set_far_nbr( far_neighbor_data *fdest,
                              int j, double d, rvec dvec )
{
  fdest->nbr = j;
  fdest->d = d;
  rvec_Copy( fdest->dvec, dvec );
  ivec_MakeZero( fdest->rel_box );
}

/* ---------------------------------------------------------------------- */

int PairReaxC::estimate_reax_lists()
{
  int itr_i, itr_j, i, j;
  int num_nbrs, num_marked;
  int *ilist, *jlist, *numneigh, **firstneigh, *marked;
  double d_sqr;
  rvec dvec;
  double **x;

  int mincap = system->mincap;
  double safezone = system->safezone;

  x = atom->x;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  num_nbrs = 0;
  num_marked = 0;
  marked = (int*) calloc( system->N, sizeof(int) );

  int numall = list->inum + list->gnum;

  for( itr_i = 0; itr_i < numall; ++itr_i ){
    i = ilist[itr_i];
    marked[i] = 1;
    ++num_marked;
    jlist = firstneigh[i];

    for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ){
      j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance( x[j], x[i], &d_sqr, &dvec );

      if (d_sqr <= SQR(control->nonb_cut))
        ++num_nbrs;
    }
  }

  free( marked );

  return static_cast<int> (MAX( num_nbrs*safezone, mincap*MIN_NBRS ));
}

/* ---------------------------------------------------------------------- */

int PairReaxC::write_reax_lists()
{
  int itr_i, itr_j, i, j;
  int num_nbrs;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double d_sqr, cutoff_sqr;
  rvec dvec;
  double *dist, **x;
  reax_list *far_nbrs;
  far_neighbor_data *far_list;

  x = atom->x;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  far_nbrs = lists + FAR_NBRS;
  far_list = far_nbrs->select.far_nbr_list;

  num_nbrs = 0;
  int inum = list->inum;
  dist = (double*) calloc( system->N, sizeof(double) );

  int numall = list->inum + list->gnum;

  for( itr_i = 0; itr_i < numall; ++itr_i ){
    i = ilist[itr_i];
    jlist = firstneigh[i];
    Set_Start_Index( i, num_nbrs, far_nbrs );

    if (i < inum)
      cutoff_sqr = control->nonb_cut*control->nonb_cut;
    else
      cutoff_sqr = control->bond_cut*control->bond_cut;

    for( itr_j = 0; itr_j < numneigh[i]; ++itr_j ){
      j = jlist[itr_j];
      j &= NEIGHMASK;
      get_distance( x[j], x[i], &d_sqr, &dvec );

      if (d_sqr <= (cutoff_sqr)) {
        dist[j] = sqrt( d_sqr );
        set_far_nbr( &far_list[num_nbrs], j, dist[j], dvec );
        ++num_nbrs;
      }
    }
    Set_End_Index( i, num_nbrs, far_nbrs );
  }

  free( dist );

  return num_nbrs;
}

/* ---------------------------------------------------------------------- */

void PairReaxC::read_reax_forces(int /*vflag*/)
{
  if(flag_bb==1){
    if(count_bb_timesteps<MAX_NUM_TIMESTEPS){
      count_bb_timesteps++;
      //add_bb_potential(); //no need
    }
    else{
      flag_bb=0;
      count_bb_timesteps=0;
      if(energy_fp!=NULL)
        fprintf(energy_fp,"\nfinish");
      if(update->ntimestep>150000)  calm_down=1000;
      //OFEK
      printf("\n\n**** finish %d timesteps at timestep %d****\n\n", MAX_NUM_TIMESTEPS, update->ntimestep); 
      }
  }

  //printf("\nadding the force");
  for( int i = 0; i < system->N; ++i ) {
    system->my_atoms[i].f[0] = workspace->f[i][0];
    system->my_atoms[i].f[1] = workspace->f[i][1];
    system->my_atoms[i].f[2] = workspace->f[i][2];
    atom->f[i][0] += -workspace->f[i][0];
    atom->f[i][1] += -workspace->f[i][1];
    atom->f[i][2] += -workspace->f[i][2];
  }

}

/* ---------------------------------------------------------------------- */

void *PairReaxC::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"chi") == 0 && chi) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) chi[i] = system->reax_param.sbp[map[i]].chi;
      else chi[i] = 0.0;
    return (void *) chi;
  }
  if (strcmp(str,"eta") == 0 && eta) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) eta[i] = system->reax_param.sbp[map[i]].eta;
      else eta[i] = 0.0;
    return (void *) eta;
  }
  if (strcmp(str,"gamma") == 0 && gamma) {
    for (int i = 1; i <= atom->ntypes; i++)
      if (map[i] >= 0) gamma[i] = system->reax_param.sbp[map[i]].gamma;
      else gamma[i] = 0.0;
    return (void *) gamma;
  }
  return NULL;
}

/* ---------------------------------------------------------------------- */

double PairReaxC::memory_usage()
{
  double bytes = 0.0;

  // From pair_reax_c
  bytes += 1.0 * system->N * sizeof(int);
  bytes += 1.0 * system->N * sizeof(double);

  // From reaxc_allocate: BO
  bytes += 1.0 * system->total_cap * sizeof(reax_atom);
  bytes += 19.0 * system->total_cap * sizeof(double);
  bytes += 3.0 * system->total_cap * sizeof(int);

  // From reaxc_lists
  bytes += 2.0 * lists->n * sizeof(int);
  bytes += lists->num_intrs * sizeof(three_body_interaction_data);
  bytes += lists->num_intrs * sizeof(bond_data);
  bytes += lists->num_intrs * sizeof(dbond_data);
  bytes += lists->num_intrs * sizeof(dDelta_data);
  bytes += lists->num_intrs * sizeof(far_neighbor_data);
  bytes += lists->num_intrs * sizeof(hbond_data);

  if(fixspecies_flag)
    bytes += 2 * nmax * MAXSPECBOND * sizeof(double);

  //mine
  bytes += 3.0 * atom->nlocal * sizeof(double); //for f_fourset
  bytes += 4.0 * atom->nlocal * sizeof(int); //for fourset

  return bytes;
}

/* ---------------------------------------------------------------------- */

void PairReaxC::FindBond()
{
  int i, j, pj, nj;
  double bo_tmp, bo_cut;

  bond_data *bo_ij;
  bo_cut = 0.10;

  for (i = 0; i < system->n; i++) {
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
        if (nj > MAXSPECBOND) error->all(FLERR,"Increase MAXSPECBOND in reaxc_defs.h");
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
//return 0 for success. else, return -1.
int PairReaxC::set_extra_potential_parameters(){
  //printf("\nPairReaxC:; in set_extra_potential_parameters\n");

//ofek
  parameters_fp = fopen("Extra_Potential_Parameters.txt","r");
  int rtn_val;

  if (parameters_fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open fix reax/c/bonds file Extra_Potential_Parameters.txt");
    error->one(FLERR,str);
    return -1;
  }
  char buff[1000];
  fread(buff, 1000, 1, parameters_fp);

  char *token = strtok(buff, "\n");
  int finish_flag=0;
  while(token){
    //printf("\ntoken %s\n", token);
    if(strcmp(token, "O-C pair tags")==0 || strcmp(token, "o-c pair tags")==0){
      for(int i=0; i<4; i++)
        token = strtok(NULL, "\n");
    }
    else if(strcmp(token, "n tags")==0 || strcmp(token, "N tags")==0){
        token = strtok(NULL, "\n");
    }
    else if(strcmp(token, "max_iterarions_of_potential")==0){
      token = strtok(NULL, "\n");
      //printf("\n max iter token= %s\n",token);
      rtn_val=sscanf(token, "%d", &MAX_NUM_TIMESTEPS);
      if(rtn_val<=0){
        fclose(parameters_fp);
        return -1;
      }
    }
    else if(strcmp(token, "TYPE1 TYPE2 F1 F2 R12")==0){
      for(int i=0; i<4; i++){
        int type1=0, type2=0;
        double temp, _f1, _f2, _r12;
        for(int j=0; j<5; j++){
          if(j<4)
            token = strtok(NULL, " ");
          else
            token = strtok(NULL, "\n");
          rtn_val=sscanf(token, "%lf", &temp);
          if(rtn_val<=0){
            fclose(parameters_fp);
            return -1;
          }
          switch(j) {
            case 0: type1=int(temp);
              break;
            case 1: type2=int(temp);
              break;
            case 2: _f1=temp;
              break;
            case 3: _f2=temp;
              break;
            case 4: _r12=temp;
              break;
          }
          if(type1==0 ||type2==0) continue;

          F1[type1][type2]=F1[type2][type1]=_f1*0.5;
          F2[type1][type2]=F2[type2][type1]=_f2;
          wanted_dist[type1][type2]=wanted_dist[type2][type1]=_r12;
          //FOR DEBUGGING
          //printf("\n~~~~type1 %d type2 %d _f1=%f _f2=%f _r12=%f MAX_NUM_TIMESTEPS=%d\n",type1, type2, _f1,_f2, _r12, MAX_NUM_TIMESTEPS);
        }
      }
      finish_flag=1;
    }
    else{
      printf("\n format error\n");
      fclose(parameters_fp);
      return -1;
    }
    if(finish_flag==1){
      fclose(parameters_fp);
      return 0;
    }
    token = strtok(NULL, "\n");
  }

  if(finish_flag==0){
    fclose(parameters_fp);
    printf("\n not finish\n");
    return -1;
  }
}

/* ---------------------------------------------------------------------- */
/*returns 1 if apply the extra potential on the foursets. else, 0*/

int PairReaxC::set_fourset(int **foursets, int num_foursets){
  if(count_bb_timesteps>0 || flag_bb==1)
    return 0;
  if(update->ntimestep<10000){
    //printf("\ntoo early\n");
    return 0;
  }
  if(calm_down>0){
    calm_down--;
    return 0;
  }

  if(update->laststep-update->ntimestep<MAX_NUM_TIMESTEPS)
  {
    //printf("\ntoo late\n");
    return 0;
  }
  int set_params_flag=set_extra_potential_parameters();
  if(set_params_flag==-1){
    printf("\nerror in set_extra_potential_parameters\n");
     error->all(FLERR,"Illegal extra_potential_parameters file reax/c command");
    return 0;
  }
  //OFEK
  printf("\n~~~in set_fourset timestep %d~~~\n", update->ntimestep);
  for(int i=0; i<num_foursets; i++){
    printf("fourset #%d: %d %d %d %d\n",i,foursets[i][0],foursets[i][1],foursets[i][2],foursets[i][3]);
  }
  if(energy_fp != NULL)
    fprintf(energy_fp,"\nstart");

  count_bb_timesteps=0;
  flag_bb=1;
  int i;
  num_fourset=num_foursets;
  //copy the f 2D array (to avoid aliasing)
  for(i=0; i<num_foursets; ++i)
    for(int j=0; j<4; j++)
      fourset[i][j]=foursets[i][j];
  for(i; i<atom->nlocal; i++)
    for(int j=0; j<4; j++)
      fourset[i][j]=0;
  return 1;
}

/* ---------------------------------------------------------------------- */
double PairReaxC::compute_BB(){
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
}

/* ---------------------------------------------------------------------- */
double PairReaxC::compute_BB_pair(int i_tag, int j_tag){
  
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
  for( pk = Start_Index(i, lists+BONDS); pk < End_Index(i, lists+BONDS); ++pk ) {
      bo_ij = &( (lists+BONDS)->select.bond_list[pk] );
      k = bo_ij->nbr;
      if (j == k || atom->tag[k] == j_tag){
        cal_dist_flag=1;
        rij=bo_ij->d;
        dvec[0]=bo_ij->dvec[0];
        dvec[1]=bo_ij->dvec[1];
        dvec[2]=bo_ij->dvec[2];
        break;
      }
  }

  if(rij==0){
    far_neighbor_data *nbr_ij;
    for( pk = Start_Index(i, lists+FAR_NBRS); pk < End_Index(i, lists+FAR_NBRS); ++pk ) {
      nbr_ij = &( (lists+FAR_NBRS)->select.far_nbr_list[pk] );
      k = nbr_ij->nbr;
      if (j == k || atom->tag[k] == j_tag){
        cal_dist_flag=2;
        rij=nbr_ij->d;
        dvec[0]=nbr_ij->dvec[0];
        dvec[1]=nbr_ij->dvec[1];
        dvec[2]=nbr_ij->dvec[2];
        break;
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
    
  /*//calculate the F force vector for atom i
  f_fourset[i_tag-1][0] += (delx*fpair)/rij;
  f_fourset[i_tag-1][1] += (dely*fpair)/rij;
  f_fourset[i_tag-1][2] += (delz*fpair)/rij;
  //calculate the F force vector for atom j       
  f_fourset[j_tag-1][0] -= (delx*fpair)/rij;
  f_fourset[j_tag-1][1] -= (dely*fpair)/rij;
  f_fourset[j_tag-1][2] -= (delz*fpair)/rij;*/

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

}
/* ---------------------------------------------------------------------- */
//this method gets:
//1. tag, type of atom i and tag, type of atom j
//2. the R(i,j)=the distance between them
//returns the calculated force.
double PairReaxC::single_BB(int i, int j, int itag, int jtag, int itype, int jtype, double rsq)
{
  double force;
  double r= rsq-wanted_dist[itype][jtype];
  double temp= -F2[itype][jtype] * r;
  force = -2 * F1[itype][jtype] * temp * exp(temp * r);
  
  //FOR DEBUGGING
  //OFEK
 // if(count_bb_timesteps==1 || MAX_NUM_TIMESTEPS-count_bb_timesteps==1){
 /*if(count_bb_timesteps%1000==0){
    printf("atomi=%d, atomj=%d, rsq=%f", i, j , rsq);
    printf("\nitype=%d, jtype=%d, F1=%f, F2=%f", itype, jtype, F1[itype][jtype], F2[itype][jtype]);
    printf("\nr=%f, temp=%f, force=%f\n\n", r, temp, force);
 }*/
  
  return force;
}

/* ---------------------------------------------------------------------- */
//NEED??
void PairReaxC::add_bb_potential(){
  printf("in add_bb_potential\n");
  if(flag_bb==0)
    return;

  for(int k = 0; k < system->N; ++k) {
    int tag=system->my_atoms[k].orig_id;
    if(tag>0){
      workspace->f[k][0]+=f_fourset[tag-1][0];
      workspace->f[k][1]+=f_fourset[tag-1][1];
      workspace->f[k][2]+=f_fourset[tag-1][2];
    }

  }
}
/* ---------------------------------------------------------------------- */

 int PairReaxC::tag_to_i(int tag){
   if( tag<0 || tag>atom->nlocal ) return -1;
   for(int i=0; i<atom->nlocal; i++){
     if(atom->tag[i] == tag) return i;
   }
   return -1;
 }