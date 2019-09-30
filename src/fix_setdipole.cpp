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

#include <string.h>
#include "fix_setdipole.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"


using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixSetDipole::FixSetDipole(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), iterstr(NULL), list(NULL)
{
  if (narg < 7 ) error->all(FLERR,"Illegal fix setdipole command");
  xstr = ystr = zstr = iterstr = NULL;
  
  // field vector component
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }
  // End of field vector components
  // max_iter
  if (strstr(arg[6],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    iterstr = new char[n];
    strcpy(iterstr,&arg[6][2]);
  } else if (strcmp(arg[6],"NULL") == 0) {
    iterstyle = NONE;
  } else {
    itervalue = force->numeric(FLERR,arg[6]);
    iterstyle = CONSTANT;
  }

  maxatom = 1;
  memory->create(sforce,maxatom,3,"setforce:sforce");
  
}

/* ---------------------------------------------------------------------- */

FixSetDipole::~FixSetDipole()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] iterstr;
  memory->destroy(sforce);

}

/* ---------------------------------------------------------------------- */

int FixSetDipole::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetDipole::init()
{
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix setdipole does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix setdipole does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix setdipole is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix setdipole does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix setdipole is invalid style");
  }
  // max_iter
  if (iterstr) {
    itervar = input->variable->find(iterstr);
    if (itervar < 0)
      error->all(FLERR,"Variable name for fix setdipole does not exist");
    if (input->variable->equalstyle(itervar)) iterstyle = EQUAL;
    else if (input->variable->atomstyle(itervar)) iterstyle = ATOM;
    else error->all(FLERR,"Variable for fix setdipole is invalid style");
  }
  
  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM || iterstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL || iterstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;
  
  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
    if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
    if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
    if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  }
  if (flag)
    error->all(FLERR,"Cannot use non-zero forces in an energy minimization");

  // Create neighbor request
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
}

/* ---------------------------------------------------------------------- */
void FixSetDipole::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- 

void FixSetDipole::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}
*/
void FixSetDipole::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ----------------------------------------------------------------------

void FixSetDipole::min_setup(int vflag)
{
  post_force(vflag);
}

---------------------------------------------------------------------- */

void FixSetDipole::post_force(int vflag)
{
  int inum, ii, jnum, jj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double **x = atom->x;
  double **mu = atom->mu;
  double *xi = atom->xi;
  double *radius = atom->radius;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double volume;
  double pi43 = 4.0/3.0*3.14159265358979323846;
  double pi4_inv = 1/4.0/3.14159265358979323846;
  double B [3];
  double r [3];
  double r1, r3, r5, mudotr;
	  
  // The neighbor list is built when it is instantiated
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"setforce:sforce");
  }

	if (varflag == CONSTANT) 
	{
	  for (int i = 0; i < nlocal; i++)
		if (mask[i] & groupbit) 
		{
/* 			printf("value = [%f,%f,%f]\n",xvalue,yvalue,zvalue);
			printf("var = [%f,%f,%f]\n",xvar,xvar,xvar);
			printf("mu = [%f,%f,%f]\n",mu[i][0],mu[i][1],mu[i][2]);

			
			if (xstyle) mu[i][0] = xvalue*volume*xi[i];
			if (ystyle) mu[i][1] = yvalue*volume*xi[i];
			if (zstyle) mu[i][2] = zvalue*volume*xi[i]; */
		}
	} else 
	{

		modify->clearstep_compute();

		if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
		else if (xstyle == ATOM)
		  input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
		if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
		else if (ystyle == ATOM)
		  input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
		if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
		else if (zstyle == ATOM)
		  input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);
        if (iterstyle == EQUAL) itervalue = input->variable->compute_equal(itervar);
		else if (iterstyle == ATOM)
		  input->variable->compute_atom(itervar,igroup,&sforce[0][2],3,0);

		modify->addstep_compute(update->ntimestep + 1);
		
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit) {

 				if (xstyle == ATOM) xvalue = sforce[i][0];
				if (ystyle == ATOM) yvalue = sforce[i][1];
				if (zstyle == ATOM) zvalue = sforce[i][2];
                if (iterstyle == ATOM) itervalue = sforce[i][2];
			}
		
	}
	
	int n_neighbors = 0;
	
	//printf("field is B = [%2.2g,%2.2g,%2.2g] \n",xvalue,yvalue,zvalue);

	for (int i=0; i<inum; i++) {
		
		if (mask[i] & groupbit) {
			
			volume = pi43*radius[i]*radius[i]*radius[i];
			
			mu[i][0] = xvalue*volume*xi[i];
			mu[i][1] = yvalue*volume*xi[i];
			mu[i][2] = zvalue*volume*xi[i];
			
			ii = ilist[i];
			jlist = firstneigh[ii];
			jnum = numneigh[ii];
			
			for (int j=0; j<jnum; j++)
			{
				jj = jlist[j];
				if (jj>n_neighbors) {n_neighbors = jj;}
			}
		}
		//if (i%12==0){
			//printf("particle %u has initial m = [%2.2g,%2.2g,%2.2g]\n",tag[i],mu[i][0],mu[i][1],mu[i][2]);
		//}
	}
	
	if(itervalue>0){
		double mu_correct [n_neighbors][3];
		//printf("There are %u particles: \n", n_neighbors);
		for (int i=0; i<n_neighbors; i++){
			
			volume = pi43*radius[i]*radius[i]*radius[i];
			
			mu_correct[i][0] = xvalue*volume*xi[i];
			mu_correct[i][1] = yvalue*volume*xi[i];
			mu_correct[i][2] = zvalue*volume*xi[i];
			//if (i%12==0){
				
				//printf("particle tag[%u]=%u has xi=%2.2f, vol = %2.2g and initially mc = [%2.2g,%2.2g,%2.2g] \n",i,tag[i],xi[i],volume,mu_correct[i][0],mu_correct[i][1],mu_correct[i][2]);
			//}
		}

		//printf("iteration...");	
		for (int iter=0; iter<itervalue; iter++){
			//printf("%u ",iter);
			for (int i=0; i<nlocal; i++){
				if (mask[i] & groupbit){
					// i is the current atom.
					// ilist is the directory of nieighbors
					
					ii = ilist[i];
					jlist = firstneigh[ii];
					jnum = numneigh[ii];
		
					B[0] = 0;
					B[1] = 0;
					B[2] = 0;
		
					// Calculate the field at this point from all other colloids
					//printf("current atom is %u ",i);
					//printf("with %u neighbors: \n",jnum);
					for (int j=0; j<jnum; j++)
					{
						if (mask[j] & groupbit) {
							jj = jlist[j];

							jj &= NEIGHMASK;
							//printf("%u:%u,%2.2f ",j,jj,mu_correct[jj][0]);
							r[0] = x[jj][0]-x[ii][0];
							r[1] = x[jj][1]-x[ii][1];
							r[2] = x[jj][2]-x[ii][2];
						
							r1 = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
							//printf("d = %2.2f, ",r1);
							r3 = r1*r1*r1;
							r5 = r3*r1*r1;
							
							mudotr = mu_correct[jj][0]*r[0] + mu_correct[jj][1]*r[1] + mu_correct[jj][2]*r[2];
							//printf("m.r = [%g,%g,%g].[%g,%g,%g]=%g \n",mu_correct[jj][0], mu_correct[jj][1], mu_correct[jj][2],r[0],r[1],r[2],mudotr);
							//printf("m.r = %f \n",mudotr);
							//if ((i==0||i==60) & j==(jnum-1)){
								//printf("\n neighbor %u:%u : m.r = [%f,%f,%f]",jj,mudotr);
								//printf("m_c = [%f,%f,%f] \n",mu_correct[i][0],mu_correct[i][1],mu_correct[i][2]);
							//}
							B[0] += pi4_inv * (3*r[0]*mudotr/r5-mu_correct[jj][0]/r3);
							B[1] += pi4_inv * (3*r[1]*mudotr/r5-mu_correct[jj][1]/r3);
							B[2] += pi4_inv * (3*r[2]*mudotr/r5-mu_correct[jj][2]/r3);
						}
					}
					//volume = pi43*radius[i]*radius[i]*radius[i];

					mu_correct[i][0] = mu[i][0]+B[0]*volume*xi[i];
					mu_correct[i][1] = mu[i][1]+B[1]*volume*xi[i];
					mu_correct[i][2] = mu[i][2]+B[2]*volume*xi[i];
					
					
					//printf("B_c = [%f,%f,%f] ",B[0],B[1],B[2]);
					//printf("m_c = [%f,%f,%f] ",mu_correct[i][0],mu_correct[i][1],mu_correct[i][2]);
					//printf("\n");
				}
			}
		}
		
		for (int i=0; i<nlocal; i++){
			//printf("setting corrected m of %u",i);
			if (mask[i] & groupbit){

				mu[i][0] = mu_correct[i][0];
				mu[i][1] = mu_correct[i][1];
				mu[i][2] = mu_correct[i][2];
			}
		}
	}
	// set dipole magnitude
	for (int i = 0; i<nlocal; i++) {
		if (mask[i] & groupbit) {
			mu[i][3] = sqrt(mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2]);
		}
	}
	//printf("\n ----- end timestep ------ \n ");
		
}

/* ---------------------------------------------------------------------- */

void FixSetDipole::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetDipole::min_post_force(int vflag)
{
  post_force(vflag);
}
