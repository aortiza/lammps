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

#ifdef PAIR_CLASS

PairStyle(biharmonic,PairBiHarmonic)

#else

#ifndef LMP_PAIR_BIHARMONIC
#define LMP_PAIR_BIHARMONIC

#include "pair.h"

namespace LAMMPS_NS {

class PairBiHarmonic : public Pair {
 public:
  PairBiHarmonic(class LAMMPS *);
  virtual ~PairBiHarmonic();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_bh_global;
  double **cut_bh,**cut_bhsq;
  double **k_outer,**k_inner;
  double **eccent;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_style command

Self-explanatory.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dipole/cut requires atom attributes q, mu, torque

The atom style defined does not have these attributes.

*/
