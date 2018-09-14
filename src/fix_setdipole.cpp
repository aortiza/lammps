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
  xstr(NULL), ystr(NULL), zstr(NULL)
{
  if (narg < 8 ) error->all(FLERR,"Illegal fix setdipole command");
  xstr = ystr = zstr = NULL;
  
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
  // cuttoff
  if (strstr(arg[6],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    cuttoffstr = new char[n];
    strcpy(cuttoffstr,&arg[6][2]);
  } else if (strcmp(arg[6],"NULL") == 0) {
    cuttoffstyle = NONE;
  } else {
    cuttoffvalue = force->numeric(FLERR,arg[6]);
    cuttoffstyle = CONSTANT;
  }
  // max_iter
  if (strstr(arg[7],"v_") == arg[7]) {
    int n = strlen(&arg[7][2]) + 1;
    iterstr = new char[n];
    strcpy(iterstr,&arg[7][2]);
  } else if (strcmp(arg[7],"NULL") == 0) {
    iterstyle = NONE;
  } else {
    itervalue = force->numeric(FLERR,arg[7]);
    itervalue = force->numeric(FLERR,arg[7]);
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
  delete [] cuttoffstr;
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
  
  // cuttoff
  if (cuttoffstr) {
    cuttoffvar = input->variable->find(cuttoffstr);
    if (cuttoffvar < 0)
      error->all(FLERR,"Variable name for fix setdipole does not exist");
    if (input->variable->equalstyle(cuttoffvar)) cuttoffstyle = EQUAL;
    else if (input->variable->atomstyle(cuttoffvar)) cuttoffstyle = ATOM;
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
  
  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
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
/*   neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1; */
  
}

/* ---------------------------------------------------------------------- */
void FixSetDipole::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

void FixSetDipole::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetDipole::post_force(int vflag)
{
  int inum, ii, jnum;
  int *ilist,*numneigh,**firstneigh;
  double **mu = atom->mu;
  double *xi = atom->xi;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double volume;
  double pi43 = 4.0/3.0*3.14159265358979323846;
  
  // The neighbor list is built when it is instantiated
  printf("creating lists ... ");
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  printf("lists created \n");

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

		modify->addstep_compute(update->ntimestep + 1);
		
		for (int i = 0; i < nlocal; i++)
			if (mask[i] & groupbit) {
/* 				printf("value = [%f,%f,%f]\n",xvalue,yvalue,zvalue);
				printf("var = [%f,%f,%f]\n",xvar,xvar,xvar);
				printf("mu = [%f,%f,%f]\n",mu[i][0],mu[i][1],mu[i][2]);
*/				
				//printf("sforce=%e ,zvalue = %e \n",sforce[i][2],zvalue);
 				if (xstyle == ATOM) xvalue = sforce[i][0];
				if (ystyle == ATOM) yvalue = sforce[i][1];
				if (zstyle == ATOM) zvalue = sforce[i][2];
			}
		
	}
	
	printf("list");
	for (int i=0; i<nlocal; i++)
		if (mask[i] & groupbit) {
			
			volume = pi43*radius[i]*radius[i]*radius[i];
			
			mu[i][0] = xvalue*volume*xi[i];
			mu[i][1] = yvalue*volume*xi[i];
			mu[i][2] = zvalue*volume*xi[i];
			
			printf("building list\n");
			
			ii = ilist[i];
			jnum = numneigh[ii];
			
			printf("inum = %u\n",inum);
			//for (int j; j<jnum; j++)
			//printf("number of neighbors = %u",jnum);
		
			mu[i][3] = sqrt(mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2]);
			printf("mu[%u] = [%f,%f,%f]\n",i,mu[i][0],mu[i][1],mu[i][2]);
		}
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
