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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "pair_biharmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBiHarmonic::PairBiHarmonic(LAMMPS *lmp) : Pair(lmp)
{
	single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairBiHarmonic::~PairBiHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_bh);
    memory->destroy(cut_bhsq);
    memory->destroy(k_outer);
    memory->destroy(k_inner);
    memory->destroy(eccent);
	
  }
}

/* ---------------------------------------------------------------------- */

void PairBiHarmonic::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,epair,fpair,ftot,fx,fy,fz;
  double rsq,r,dr,rk;
  double upllx, uplly, upllz, uprpx, uprpy, uprpz;
  double delpll, delprp, delprpsq, fpll, fprp, sign_delpll, abs_delpll, Epll;
	  
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  ev_init(eflag,vflag);
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
    
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
		
      j = jlist[jj];
	  //factor_lj = special_lj[sbmask(j)];
      //factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;
      delx = - xtmp + x[j][0];
      dely = - ytmp + x[j][1];
      delz = - ztmp + x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
	  
	  // Vector of direction of trap
	  dr = sqrt(mu[j][0]*mu[j][0] + mu[j][1]*mu[j][1] + mu[j][2]*mu[j][2]);
	  // Unit vector of direction of the trap 
	  if (dr>0){
		  upllx = mu[j][0]/dr;
		  uplly = mu[j][1]/dr;
		  upllz = mu[j][2]/dr;
	  } 
	  else {
		  upllx = 1;
		  uplly = 0;
		  upllz = 0;
	  }
	  // Unit vector of direction perpendicular to the trap (radial cylindrical)
	  // The case in which the particle sits exactly in the vector del, is a 
	  // fringe case and should be treated separately to assign an arbitrary 
	  // perpendicular direction
	  
	  // component of position parallel to the trap
	  delpll = delx*upllx + dely*uplly + delz*upllz;
	  // Component of position perpendicular to the trap
	  delprpsq = rsq-delpll*delpll;
	  delprp = sqrt(delprpsq);
	  
	  if (delprp>0){
		  uprpx = (delx - delpll*upllx)/delprp;
		  uprpy = (dely - delpll*uplly)/delprp;
		  uprpz = (delz - delpll*upllz)/delprp;
	  }
	  else {
		  uprpx = 0;
		  uprpy = 1;
		  uprpz = 0;
	  }

	  dr = dr/2;
	  
	  /*printf("cut = %f, d1 = %f, d2 = %f, d3 = %f, is?%u \n",
	  	cut_bh[itype][jtype],
		sqrt(pow(delpll-dr,2)+delprpsq),
		sqrt(pow(delpll+dr,2)+delprpsq),
		delprp,
		( (sqrt(pow(delpll-dr,2)+delprpsq) < cut_bh[itype][jtype] |
				   sqrt(pow(delpll+dr,2)+delprpsq) < cut_bh[itype][jtype]) & 
			       (delprp < cut_bh[itype][jtype])));*/

	  abs_delpll = sqrt(delpll*delpll);
	  
	  // This cutoff needs to be more sophisticated
	  /*printf("cut = %f, delprp = %f, dellob = %f, inside? = %u\n",
	      cut_bh[itype][jtype],delprp,sqrt(pow(abs_delpll-dr,2)+delprpsq),
		  ((sqrt(pow(abs_delpll-dr,2)+delprpsq) < cut_bh[itype][jtype]) & 
		  	       (delprp < cut_bh[itype][jtype])));*/
		
      if ((sqrt(pow(abs_delpll-dr,2)+delprpsq) < cut_bh[itype][jtype]) & 
	       (delprp < cut_bh[itype][jtype])) {
		
        fx = fy = fz = 0.0;
	
		sign_delpll = copysign(1,delpll);
		
		if (abs_delpll>=dr){
			fpll = -k_outer[itype][jtype]*(abs_delpll-dr)*sign_delpll;
			Epll = k_outer[itype][jtype]/2*(abs_delpll-dr)*(abs_delpll-dr);
	//		printf("Push In ");
	    }
		else if (abs_delpll<dr){
			fpll = k_inner[itype][jtype]*delpll;
			Epll = -k_inner[itype][jtype]/2*(delpll*delpll+dr*dr);
	//		printf("Push Out ");
		}
		fprp = -k_outer[itype][jtype]*delprp*eccent[itype][jtype];
			
		fx = fpll * upllx + fprp*uprpx;
		fy = fpll * uplly + fprp*uprpy;
		fz = fpll * upllz + fprp*uprpz;
		
        // force & torque accumulation

        f[i][0] -= fx;
        f[i][1] -= fy;
        f[i][2] -= fz;
		
		ftot = fx*fx+fy*fy+fz*fz;

	    // force & energy

	    if (r > 0.0) fpair = -2.0*ftot/r;
	    else fpair = 0.0;
		
	    if (eflag) epair = k_outer[itype][jtype]*delprp*delprp+Epll;
		
        if (newton_pair || j < nlocal) {
          f[j][0] += fx;
          f[j][1] += fy;
          f[j][2] += fz;
        }
		
        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 epair,0.0,fpair,delx,dely,delz);
								 
		
      }
    }
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
  
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBiHarmonic::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  
  memory->create(cut_bh,n+1,n+1,"pair:cut_lj");
  memory->create(cut_bhsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(k_outer,n+1,n+1,"pair:epsilon");
  memory->create(k_inner,n+1,n+1,"pair:sigma");
  memory->create(eccent,n+1,n+1,"pair:sigma");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBiHarmonic::settings(int narg, char **arg)
{
  if (narg == 0 || narg > 1)
    error->all(FLERR,"Incorrect args in pair_style command");

  cut_bh_global = force->numeric(FLERR,arg[0]);
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_bh[i][j] = cut_bh_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBiHarmonic::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);
  
  double eccent_one = 1;
  double k_outer_one = force->numeric(FLERR,arg[2]);
  double k_inner_one = force->numeric(FLERR,arg[3]);

  double cut_bh_one = cut_bh_global;
  if (narg >= 5) cut_bh_one = force->numeric(FLERR,arg[4]);
  //printf("k_in = %f, k_out = %f\n",k_inner_one,k_outer_one);
  
  if (narg >= 6) eccent_one = force->numeric(FLERR,arg[5]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      k_outer[i][j] = k_outer_one;
      k_inner[i][j] = k_inner_one;
      cut_bh[i][j] = cut_bh_one;
	  eccent[i][j] = eccent_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBiHarmonic::init_style()
{
  if (!atom->mu_flag)
    error->all(FLERR,"Pair biharmonic requires atom attributes mu");
  neighbor->request(this,instance_me);
  
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBiHarmonic::init_one(int i, int j)
{
    if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

    return cut_bh[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBiHarmonic::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&k_outer[i][j],sizeof(double),1,fp);
        fwrite(&k_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut_bh[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBiHarmonic::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&k_outer[i][j],sizeof(double),1,fp);
          fread(&k_inner[i][j],sizeof(double),1,fp);
          fread(&cut_bh[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&k_outer[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&k_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_bh[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBiHarmonic::write_restart_settings(FILE *fp)
{
  fwrite(&cut_bh_global,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBiHarmonic::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_bh_global,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_bh_global,1,MPI_DOUBLE,0,world);
}

double PairBiHarmonic::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj,
                         double &fforce)
{
  fforce = -k_outer[itype][jtype]*rsq;

  return fforce;
}