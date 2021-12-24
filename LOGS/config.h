//-----------------------------------------------------------
// configuration file to specify
// 2-D lattice size and 2-D parallelism
//-----------------------------------------------------------

#ifndef _CONFIG_H_
#define _CONFIG_H_

// MPI grid size:
static const int _Px=48;
static const int _Py=4;

extern int Px;
extern int Py;
static int _Pnum=_Px*_Py;

// local lattice size
extern int NNx;
extern int NNy;
static const int Nx=60;
static const int Ny=60;


// number of iteratoin for the solever --- it is a fixed iter. solver 
static const int IterForSolver=1000;

// frequency to call norm2
static const int norm2_freq=10;

// dirctions:  +x, -x, +y, -y, number of direcitons
enum Dir_pm {DirXp, DirXm, DirYp, DirYm, NDir_pm};

// additional constant in the "Laplace" eq. --- non-zero value gives actually a Helmholtz eq.
static const double mass2=0.1;

// max buffer size to be sent/recieved
//   buffer size = Factor * (needed size)  (+ size for flag)
//   Factor = 1,2,4,...,2^{FactIterMax-1}
static const int FactIterMax=14;


#endif
