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

#include <math.h>
#include <stdlib.h>
#include "bond_biharmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBiHarmonic::BondBiHarmonic(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondBiHarmonic::~BondBiHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(h);
    memory->destroy(eta);
  }
}

/* ---------------------------------------------------------------------- */

void BondBiHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,mx,my,mz,fx,fy,fz,ftot,ebond,fbond;
  double rsq,r,dr,rk;
  double upllx, uplly, upllz, uprpx, uprpy, uprpz;
  double delpll, delprp, fpll, fprp, sign_delpll, abs_delpll, Epll, hvalue, kt;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **mu = atom->mu; // This is the magnetic moment. 
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
	
	// Vector of direction of trap
	dr = sqrt(mu[i1][0]*mu[i1][0] + mu[i1][1]*mu[i1][1] + mu[i1][2]*mu[i1][2]);
	// Unit vector of direction of the trap 
	if (dr>0){
		upllx = mu[i1][0]/dr;
		uplly = mu[i1][1]/dr;
		upllz = mu[i1][2]/dr;
	} 
	else {
		upllx = 1;
		uplly = 0;
		upllz = 0;
	}
	//printf("%f, %f, %f, dr=%f\n",upllx,uplly,upllz,dr);
	dr = dr/2;

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    
	// component of position parallel to the trap
	delpll = delx*upllx + dely*uplly + delz*upllz;
	// Component of position perpendicular to the trap
	delprp = sqrt(rsq-delpll*delpll);

	// Unit vector of direction perpendicular to the trap (radial cylindrical)
	// The case in which the particle sits exactly in the vector del, is a 
	// fringe case and should be treated separately to assign an arbitrary 
	// perpendicular direction
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
	
	sign_delpll = copysign(1,delpll);
	abs_delpll = sqrt(delpll*delpll);
	
	hvalue = h[type];
	kt = k[type]*eta[type];
//	printf("abs_delpll = %f, dr = %f ", abs_delpll, dr);
//	printf("delx = %f, dely = %f , delz = %f ", delx, dely, delz);
//	printf("delprp = %f, delpll = %f , r = %f, rsq = %f", delprp, delpll,r,rsq);
	// Calculate Parallel Force and Parallel energy
	if (abs_delpll>=dr){
		fpll = -k[type]*(abs_delpll-dr)*sign_delpll;
		Epll = k[type]/2*(abs_delpll-dr)*(abs_delpll-dr);
//		printf("Push In ");
    } 
	else if (abs_delpll<dr){
		fpll = hvalue*delpll;
		Epll = -hvalue/2*(delpll*delpll+dr*dr);
//		printf("Push Out ");
	}

	fprp = -delprp*kt;
//	printf("fpll = %f, fprp = %f \n",fpll,fprp);
		
	fx = fpll * upllx + fprp*uprpx;
	fy = fpll * uplly + fprp*uprpy;
	fz = fpll * upllz + fprp*uprpz;
	
	ftot = fx*fx+fy*fy+fz*fz;

    // force & energy

    if (r > 0.0) fbond = -2.0*ftot/r;
    else fbond = 0.0;

    if (eflag) ebond = kt*delprp*delprp+Epll;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fx;
      f[i1][1] += fy;
      f[i1][2] += fz;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= fx;
      f[i2][1] -= fy;
      f[i2][2] -= fz;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondBiHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(h,n+1,"bond:h");
  memory->create(eta,n+1,"bond:eta");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBiHarmonic::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double h_one = force->numeric(FLERR,arg[2]);
  double eta_one = 1;
  
  if (narg < 3) eta_one = force->numeric(FLERR,arg[3]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    h[i] = h_one;
	eta[i] = eta_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}


/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondBiHarmonic::equilibrium_distance(int i)
{
  return 0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBiHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&h[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&eta[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBiHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&h[1],sizeof(double),atom->nbondtypes,fp);
    fread(&eta[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&h[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&eta[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondBiHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],h[i],eta[i]);
}

/* ---------------------------------------------------------------------- */

double BondBiHarmonic::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r;
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0*rk/r;
  return rk*dr;
}
