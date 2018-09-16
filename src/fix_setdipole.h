/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(setdipole,FixSetDipole)

#else

#ifndef LMP_FIX_SETDIPOLE_H
#define LMP_FIX_SETDIPOLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSetDipole : public Fix {
 public:
  FixSetDipole(class LAMMPS *, int, char **);
  ~FixSetDipole();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  //void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  
 private:
  int ilevel_respa;
	  
 protected:
  double xvalue,yvalue,zvalue;
  int itervalue;
  int varflag;
  char *xstr,*ystr,*zstr,*cuttoffstr,*iterstr;
  int xvar,yvar,zvar,cuttoffvar,itervar,xstyle,ystyle,zstyle,cuttoffstyle,iterstyle;
  
  int maxatom;
  double **sforce;

  class NeighList *list;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

Self-explanatory.

*/
